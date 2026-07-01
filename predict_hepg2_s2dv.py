import os
import pickle as pkl

# Workaround for loading S2DV pickle files when XGBoost version is incompatible.
# This is safe only if we do NOT use the XGBoost model afterwards.
try:
    import xgboost
    from xgboost.core import Booster

    def _skip_xgb_booster_setstate(self, state):
        self.handle = None
        self.__dict__.update({"_xgb_booster_skipped": True})

    Booster.__setstate__ = _skip_xgb_booster_setstate
    print("Patched XGBoost Booster.__setstate__; XGBoost models will be skipped.")
except Exception as e:
    print(f"XGBoost patch skipped: {e}")

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def mol2alt_sentence(mol, radius):
    radii = list(range(int(radius) + 1))
    info = {}
    _ = AllChem.GetMorganFingerprint(mol, radius, bitInfo=info)

    mol_atoms = [a.GetIdx() for a in mol.GetAtoms()]
    dict_atoms = {x: {r: None for r in radii} for x in mol_atoms}

    for element in info:
        for atom_idx, radius_at in info[element]:
            dict_atoms[atom_idx][radius_at] = element

    identifiers_alt = []
    for atom in dict_atoms:
        for r in radii:
            identifiers_alt.append(dict_atoms[atom][r])

    alternating_sentence = map(str, [x for x in identifiers_alt if x])
    return list(alternating_sentence)


def get_ECFP(mol, radius):
    ecfps = mol2alt_sentence(mol, radius)
    if len(ecfps) == 0:
        return []

    if len(ecfps) % (radius + 1) != 0:
        ecfps = ecfps[:-(len(ecfps) % (radius + 1))]

    if len(ecfps) == 0:
        return []

    ecfp_by_radius = list(
        np.array(ecfps).reshape(
            int(len(ecfps) / (radius + 1)),
            radius + 1
        )[:, radius]
    )
    return ecfp_by_radius


def get_sentence_vec(tokens, embedding, token_dict):
    feature_vec = np.zeros(512)
    n_missing = 0

    for token in tokens:
        if token in token_dict:
            feature_vec = np.add(feature_vec, embedding[token_dict[token]])
        else:
            n_missing += 1

    n_valid = len(tokens) - n_missing

    if n_valid <= 0:
        return np.zeros(512)

    return np.divide(feature_vec, n_valid)


def predict_one_smiles(smiles, model, token_dict, embedding, model_name="SVM"):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return {
            "smiles": smiles,
            "valid_smiles": False,
            "hepg2_pred_label": np.nan,
            "hepg2_toxic_proba": np.nan
        }

    tokens = get_ECFP(mol, 1)
    vec = get_sentence_vec(tokens, embedding, token_dict)

    selected_model = None
    for name, m in model:
        if name == model_name:
            selected_model = m
            break

    if selected_model is None:
        raise ValueError(f"Cannot find model named {model_name}. Available: {[x[0] for x in model]}")

    pred_label = selected_model.predict(vec.reshape(1, -1))[0]

    if hasattr(selected_model, "predict_proba"):
        pred_proba = selected_model.predict_proba(vec.reshape(1, -1))[:, 1][0]
    else:
        pred_proba = np.nan

    return {
        "smiles": smiles,
        "valid_smiles": True,
        "hepg2_pred_label": int(pred_label),
        "hepg2_toxic_proba": float(pred_proba)
    }


def main(input_csv, output_csv, smiles_col="smiles"):
    model_root = "./model"

    hepg2_model = pkl.load(open(os.path.join(model_root, "HepG2.ECFP.models.pkl"), "rb"))
    hepg2_token = pkl.load(open(os.path.join(model_root, "HepG2_token.pkl"), "rb"))
    hepg2_emb = pkl.load(open(os.path.join(model_root, "HepG2_emb.pkl"), "rb"))

    # df = pd.read_csv(input_csv)
    df = pd.read_excel(input_csv)


    if smiles_col not in df.columns:
        raise ValueError(f"Cannot find SMILES column: {smiles_col}. Existing columns: {list(df.columns)}")

    results = []
    for smi in df[smiles_col].astype(str).tolist():
        res = predict_one_smiles(
            smiles=smi,
            model=hepg2_model,
            token_dict=hepg2_token,
            embedding=hepg2_emb,
            model_name="SVM"
        )
        results.append(res)

    pred_df = pd.DataFrame(results)
    out = pd.concat([df.reset_index(drop=True), pred_df.drop(columns=["smiles"])], axis=1)
    out.to_csv(output_csv, index=False)
    print(f"Saved predictions to: {output_csv}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--smiles_col", default="smiles")
    args = parser.parse_args()

    main(args.input_csv, args.output_csv, args.smiles_col)