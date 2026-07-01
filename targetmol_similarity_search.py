# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd
from tqdm import tqdm

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize


def read_table(path):
    ext = os.path.splitext(path)[1].lower()

    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path)

    if ext == ".csv":
        for enc in ["utf-8-sig", "utf-8", "gbk", "gb18030"]:
            try:
                return pd.read_csv(path, encoding=enc)
            except Exception:
                pass
        raise ValueError("CSV 读取失败，请检查编码。")

    raise ValueError(f"不支持的文件格式：{ext}")


def standardize_smiles(smiles):
    """
    SMILES 标准化：
    1. RDKit 解析
    2. Cleanup
    3. 取最大分子片段，常用于去盐
    4. 去电荷
    5. 输出 canonical SMILES
    """
    if pd.isna(smiles):
        return None

    smiles = str(smiles).strip()
    if smiles == "" or smiles.lower() in ["nan", "none", "null"]:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = rdMolStandardize.Cleanup(mol)

        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)

        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)

        Chem.SanitizeMol(mol)

        return Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=True
        )

    except Exception:
        return None


def smiles_to_fp(smiles, n_bits=2048):
    if smiles is None or pd.isna(smiles):
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # ECFP4: Morgan fingerprint radius = 2
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius=2,
        nBits=n_bits
    )

    return fp


def prepare_seed_df(seed_file, smiles_col="SMILES"):
    seed_df = read_table(seed_file)

    if smiles_col not in seed_df.columns:
        raise ValueError(f"种子文件中没有找到 SMILES 列：{smiles_col}")

    seed_df = seed_df[[smiles_col]].copy()
    seed_df = seed_df.rename(columns={smiles_col: "seed_raw_smiles"})

    seed_df["seed_std_smiles"] = [
        standardize_smiles(smi)
        for smi in tqdm(seed_df["seed_raw_smiles"], desc="Standardizing seed SMILES")
    ]

    seed_df["seed_fp"] = seed_df["seed_std_smiles"].apply(smiles_to_fp)

    before = len(seed_df)
    seed_df = seed_df[seed_df["seed_fp"].notna()].copy()
    after = len(seed_df)

    print(f"有效种子化合物数量：{after} / {before}")

    seed_df = seed_df.drop_duplicates(subset=["seed_std_smiles"]).reset_index(drop=True)

    return seed_df


def prepare_targetmol_df(lib_file, id_col="ID", smiles_col="SMILES"):
    lib_df = read_table(lib_file)

    for col in [id_col, smiles_col]:
        if col not in lib_df.columns:
            raise ValueError(f"TargetMol 文件中没有找到列：{col}")

    lib_df = lib_df[[id_col, smiles_col]].copy()
    lib_df = lib_df.rename(columns={
        id_col: "targetmol_id",
        smiles_col: "targetmol_raw_smiles"
    })

    lib_df["targetmol_std_smiles"] = [
        standardize_smiles(smi)
        for smi in tqdm(lib_df["targetmol_raw_smiles"], desc="Standardizing TargetMol SMILES")
    ]

    lib_df["targetmol_fp"] = lib_df["targetmol_std_smiles"].apply(smiles_to_fp)

    before = len(lib_df)
    lib_df = lib_df[lib_df["targetmol_fp"].notna()].copy()
    after = len(lib_df)

    print(f"有效 TargetMol 分子数量：{after} / {before}")

    lib_df = lib_df.drop_duplicates(subset=["targetmol_std_smiles"]).reset_index(drop=True)

    return lib_df


def search_similarity(seed_df, lib_df, top_k=5, threshold=0.75):
    results = []

    lib_fps = list(lib_df["targetmol_fp"])

    for seed_idx, seed_row in tqdm(
        seed_df.iterrows(),
        total=len(seed_df),
        desc="Searching similar compounds"
    ):
        seed_fp = seed_row["seed_fp"]

        similarities = DataStructs.BulkTanimotoSimilarity(seed_fp, lib_fps)

        sim_df = pd.DataFrame({
            "targetmol_index": range(len(lib_df)),
            "similarity": similarities
        })

        sim_df = sim_df[sim_df["similarity"] >= threshold]
        sim_df = sim_df.sort_values("similarity", ascending=False).head(top_k)

        for rank, (_, hit) in enumerate(sim_df.iterrows(), start=1):
            lib_row = lib_df.iloc[int(hit["targetmol_index"])]

            results.append({
                "seed_index": seed_idx,
                "seed_raw_smiles": seed_row["seed_raw_smiles"],
                "seed_std_smiles": seed_row["seed_std_smiles"],

                "targetmol_id": lib_row["targetmol_id"],
                "targetmol_raw_smiles": lib_row["targetmol_raw_smiles"],
                "targetmol_std_smiles": lib_row["targetmol_std_smiles"],

                "similarity": hit["similarity"],
                "rank_for_seed": rank
            })

    result_df = pd.DataFrame(results)

    return result_df


def deduplicate_hits(result_df):
    """
    一个 TargetMol 分子可能被多个 seed 命中。
    这里按 targetmol_id 去重，只保留 similarity 最高的一条。
    """
    if result_df.empty:
        return result_df

    dedup_df = (
        result_df
        .sort_values("similarity", ascending=False)
        .drop_duplicates(subset=["targetmol_id"], keep="first")
        .reset_index(drop=True)
    )

    dedup_df["global_rank_by_max_similarity"] = range(1, len(dedup_df) + 1)

    return dedup_df


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seed_file", required=True)
    parser.add_argument("--lib_file", required=True)

    parser.add_argument("--seed_smiles_col", default="SMILES")
    parser.add_argument("--lib_id_col", default="ID")
    parser.add_argument("--lib_smiles_col", default="SMILES")

    parser.add_argument("--top_k", type=int, default=5)
    parser.add_argument("--threshold", type=float, default=0.75)

    parser.add_argument("--output_prefix", default="TargetMol_similarity")

    args = parser.parse_args()

    seed_df = prepare_seed_df(
        seed_file=args.seed_file,
        smiles_col=args.seed_smiles_col
    )

    lib_df = prepare_targetmol_df(
        lib_file=args.lib_file,
        id_col=args.lib_id_col,
        smiles_col=args.lib_smiles_col
    )

    result_df = search_similarity(
        seed_df=seed_df,
        lib_df=lib_df,
        top_k=args.top_k,
        threshold=args.threshold
    )

    dedup_df = deduplicate_hits(result_df)

    out1 = args.output_prefix + "_per_seed_top_hits.xlsx"
    out2 = args.output_prefix + "_dedup_max_similarity.xlsx"

    result_df.to_excel(out1, index=False)
    dedup_df.to_excel(out2, index=False)

    print("Done.")
    print(f"每个 seed 的 Top hits：{out1}")
    print(f"按 TargetMol ID 去重后的候选：{out2}")
    print(f"per-seed hits 数量：{len(result_df)}")
    print(f"dedup hits 数量：{len(dedup_df)}")


if __name__ == "__main__":
    main()