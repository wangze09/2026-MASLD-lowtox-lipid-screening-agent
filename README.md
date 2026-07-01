[README.md](https://github.com/user-attachments/files/29531435/README.md)
# AI-Assisted Lipid-Lowering & Low-Toxicity Compound Screening

> **Competition:** AI for Science Competition — Life Science Track  
> **Task:** Nominate candidate molecules with lipid-lowering activity and low hepatocytotoxicity from a commercial library (TargetMol ~23k compounds) in the HepG2-FFA model

---

## Overview

This repository documents a multi-step AI-assisted virtual screening workflow designed to identify novel, low-toxic lipid-lowering compounds for the HepG2 free-fatty-acid (FFA) intracellular lipid accumulation model. The pipeline integrates literature-based seed curation, molecular fingerprint similarity search, dual-model cytotoxicity prediction, and LLM-assisted multi-criteria ranking.

---

## Workflow Summary

```
Step 1: GPT-assisted seed compound curation (literature mining)
           ↓
Step 2: ECFP4-based similarity search against TargetMol library (Tanimoto ≥ 0.75, Top-5)
           ↓
Step 3a: S2DV model — HepG2 cytotoxicity prediction (local)
Step 3b: Nature/chemprop model — HepG2 cytotoxicity prediction (server)
           ↓
Step 4: Manual integration → final_candidates_info.csv
           ↓
Step 5: LLM multi-criteria scoring & nomination → Top-10 candidates
```

---

## Repository Structure

```
.
├── targetmol_similarity_search.py       # Step 2: ECFP4 similarity search
├── S2DV-main/S2DV-main/
│   ├── predict_hepg2_s2dv.py            # Step 3a: S2DV cytotoxicity prediction wrapper
│   ├── S2DV_main.py                     # Core S2DV model code
│   └── model/                           # Pre-trained S2DV model files
├── antibioticsai/                       # Step 3b: chemprop Nature model
│   └── final_checkpoints/cytotox_hepg2/ # Pre-trained chemprop checkpoints
│
├── GPT_总结的降脂低毒种子化合物表格_wz_jw校验.xlsx   # Curated seed compounds (SMILES + evidence)
├── HepG2_FFA_TargetMol_similarity_075_dedup_max_similarity.xlsx  # Similarity search output
├── s2dv_hepg2_pred_candidates.csv       # S2DV prediction results
├── candidate_targetmol_pred_nature_HepG2.csv  # Nature chemprop prediction results
├── final_candidates_info.csv            # Merged candidate table for LLM ranking
│
├── 候选分子提名清单_Top10.csv            # Final Top-10 nomination list
├── 机制验证方案_Top10.csv               # Mechanistic hypotheses & validation plans
└── 机制与验证方案_Top10.pdf             # Submission-ready PDF report
```

---

## Step-by-Step Reproduction

### Environment

```bash
conda create -n lipid_screen python=3.9
conda activate lipid_screen
pip install -r requirements.txt
```

For the Nature/chemprop model (Step 3b), a separate server environment is required (see `antibioticsai/` for `environment.yml`).

---

### Step 1 — Seed Compound Curation

Seed compounds were identified by prompting GPT with:

> "Search the literature for candidate molecules that can reduce lipid accumulation in HepG2 cells under FFA-induced conditions, while not significantly impairing cell viability at effective concentrations. Provide SMILES or compound names."

Results were manually reviewed and cross-validated. The final curated seed list with SMILES is in:
`GPT_总结的降脂低毒种子化合物表格_wz_jw校验.xlsx`

---

### Step 2 — Similarity Search

```bash
python targetmol_similarity_search.py \
  --seed_file "GPT_总结的降脂低毒种子化合物表格_wz_jw校验.xlsx" \
  --lib_file "T001 TargetMol现货产品22966.csv" \
  --seed_smiles_col "SMILES" \
  --threshold 0.75 \
  --top_k 5 \
  --output_prefix "HepG2_FFA_TargetMol_similarity_075"
```

**Output:**
- `HepG2_FFA_TargetMol_similarity_075_per_seed_top_hits.xlsx` — all hits per seed
- `HepG2_FFA_TargetMol_similarity_075_dedup_max_similarity.xlsx` — deduplicated by max similarity

**Parameters:**
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Fingerprint | ECFP4 (Morgan r=2, 2048 bits) | Standard scaffold-similarity metric |
| Tanimoto threshold | 0.75 | Balances novelty and structural relevance |
| Top-K per seed | 5 | Limits output size while ensuring coverage |

---

### Step 3a — S2DV HepG2 Cytotoxicity Prediction

```bash
cd S2DV-main/S2DV-main
python predict_hepg2_s2dv.py \
  --input_csv ../../HepG2_FFA_TargetMol_similarity_075_dedup_max_similarity.xlsx \
  --output_csv ../../s2dv_hepg2_pred_candidates.csv \
  --smiles_col targetmol_std_smiles
```

**Model source:** [NTU-MedAI/S2DV](https://github.com/NTU-MedAI/S2DV)  
**Output column:** `hepg2_toxic_proba` (SVM classifier, probability of HepG2 toxicity)

---

### Step 3b — Nature/chemprop HepG2 Cytotoxicity Prediction

Run on server with chemprop environment:

```bash
chemprop_predict \
  --test_path candidate_targetmol.csv \
  --checkpoint_dir antibioticsai/final_checkpoints/cytotox_hepg2 \
  --preds_path candidate_targetmol_pred_nature_HepG2.csv \
  --features_generator rdkit_2d_normalized \
  --no_features_scaling
```

**Model source:** [felixjwong/antibioticsai](https://github.com/felixjwong/antibioticsai)  
**Output column:** `TOXICITY` (probability of HepG2 toxicity)

---

### Step 4 — Manual Data Integration

The two prediction outputs were merged with similarity search results and TargetMol compound annotations to produce `final_candidates_info.csv`.

**Key columns in `final_candidates_info.csv`:**

| Column | Description |
|--------|-------------|
| `seed_name` | Name of the matched seed compound |
| `HepG2-FFA模型证据` | Literature evidence for lipid-lowering in HepG2-FFA model |
| `targetmol_id` | TargetMol compound ID for submission |
| `similarity` | ECFP4 Tanimoto similarity to seed |
| `hepg2_toxic_proba_S2DV` | HepG2 toxicity probability (S2DV model) |
| `hepg2_toxic_proba_nature` | HepG2 toxicity probability (Nature chemprop model) |
| `Name`, `Pathways`, `Target`, `Bioactivity` | Vendor-provided annotations |

---

### Step 5 — LLM-Assisted Multi-Criteria Ranking

`final_candidates_info.csv` was submitted to an LLM with the following scoring criteria:

- **Lipid-lowering potential** — seed evidence + structural features + pathway annotations
- **Toxicity assessment** — mean of both predicted probabilities, threshold < 0.35 preferred
- **Novelty** — Tanimoto similarity to seed (lower = more novel)
- **Mechanistic plausibility** — SREBP1c/FASN/ACC (de novo lipogenesis), PPARα/AMPK/CPT1 (β-oxidation), autophagy

The strategy was **5 confident hits (similarity ≥ 0.9) + 5 novel exploratory candidates (similarity 0.75–0.85)**.

**Outputs:**
- `候选分子提名清单_Top10.csv` — ranked Top-10 with scores and rationale
- `机制验证方案_Top10.csv` — mechanistic hypotheses and validation plans
- `机制与验证方案_Top10.pdf` — PDF submission

---

## Key Results

| Rank | Compound | TargetMol ID | Strategy | Similarity | Mean Tox Prob | Score |
|------|----------|--------------|----------|------------|---------------|-------|
| 1 | Kaempferol | T2177 | Confident | 1.00 | 0.191 | 85.2 |
| 2 | p-Coumaric acid | T2863 | Confident | 1.00 | 0.181 | 83.5 |
| 3 | Morin | T2835 | Exploratory | 0.75 | 0.201 | 82.0 |
| 4 | Curcumin | T1516 | Confident | 1.00 | 0.205 | 81.9 |
| ... | ... | ... | ... | ... | ... | ... |

---

## External Models & References

| Model | Repository | Paper |
|-------|------------|-------|
| S2DV (HepG2 cytotoxicity) | [NTU-MedAI/S2DV](https://github.com/NTU-MedAI/S2DV) | — |
| Chemprop cytotoxicity | [felixjwong/antibioticsai](https://github.com/felixjwong/antibioticsai) | Wong et al., *Nature* |
| Compound library | TargetMol T001 (~22,966 compounds) | — |

---

## License

This repository is for competition submission and academic reference only.
The pre-trained model weights in `S2DV-main/` and `antibioticsai/` are subject to their respective upstream licenses.
