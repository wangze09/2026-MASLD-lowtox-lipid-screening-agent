# AI-Assisted Lipid-Lowering & Low-Toxicity Compound Screening

> **Competition:** AI for Science Competition — Life Science Track  
> **Task:** Nominate candidate molecules with lipid-lowering activity and low hepatocytotoxicity from a commercial library (TargetMol ~23k compounds) in the HepG2-FFA model

---

## Overview

This repository provides a reproducible AI-assisted screening workflow for the MASLD low-toxicity lipid-lowering compound discovery task. The workflow is organized as a modular screening agent rather than a single predictive model. It integrates literature-based seed compound curation, ECFP4/Tanimoto similarity expansion against the TargetMol compound library, dual HepG2 cytotoxicity prediction, manual evidence integration, and LLM-assisted multi-criteria ranking.

The final nomination strategy prioritizes compounds that balance lipid-lowering potential, predicted low hepatocytotoxicity, mechanistic plausibility, and structural novelty. The submitted Top-10 list contains both high-confidence candidates with strong seed-compound support and exploratory candidates with lower structural similarity but plausible lipid-metabolism mechanisms.

---

## Workflow Summary

```
Step 1: LLM-assisted seed compound curation (literature mining)
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
├── competition_library_structures.csv   # Structure library for Step 2 (ID/CAS/SMILES/Formula/MolWt)
├── predict_hepg2_s2dv.py            # Step 3a: S2DV cytotoxicity prediction wrapper
├── External dependency: felixjwong/antibioticsai                      # Step 3b: chemprop Nature model
│   └── final_checkpoints/cytotox_hepg2/ # Pre-trained chemprop checkpoints(Reproduce following to the original Github)
├── Agent总结的降脂低毒种子化合物表格_校验版.xlsx   # Curated seed compounds (SMILES + evidence)
├── HepG2_FFA_TargetMol_similarity_075_dedup_max_similarity.xlsx  # Similarity search output
├── s2dv_hepg2_pred_candidates.csv       # S2DV prediction results
├── candidate_targetmol_pred_nature_HepG2.csv  # Nature chemprop prediction results
├── final_candidates_info.csv            # Merged candidate table for LLM ranking
│
├── 候选分子提名清单_Top10.csv            # Final Top-10 nomination list    (Upload in official channel)
├── 机制验证方案_Top10.csv               # Mechanistic hypotheses & validation plans (Upload in official channel)
└── 机制与验证方案_Top10.pdf             # Submission-ready PDF report  (Upload in official channel)
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
`Agent总结的降脂低毒种子化合物表格_校验版.xlsx`

---

### Step 2 — Similarity Search

```bash
python targetmol_similarity_search.py \
  --seed_file "Agent总结的降脂低毒种子化合物表格_校验版.xlsx" \
  --lib_file "competition_library_structures.csv" \
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
python predict_hepg2_s2dv.py \
  --input_csv ../../HepG2_FFA_TargetMol_similarity_075_dedup_max_similarity.xlsx \
  --output_csv ../../s2dv_hepg2_pred_candidates.csv \
  --smiles_col targetmol_std_smiles
```

**Model source:** [NTU-MedAI/S2DV](https://github.com/NTU-MedAI/S2DV)  
**Output column:** `hepg2_toxic_proba` (SVM classifier, probability of HepG2 toxicity)

### S2DV pretrained model files

The pretrained S2DV model files are not redistributed in this repository.
Please obtain the following files from the official NTU-MedAI/S2DV repository
and place them in the `model/` directory:

- HepG2.ECFP.models.pkl
- HepG2_token.pkl
- HepG2_emb.pkl
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

The two prediction outputs were merged with similarity search results to produce `final_candidates_info.csv`. For candidate-level annotation, an agent was used to query the corresponding compounds on the TargetMol website by CAS number or SMILES, and to collect the matched compound name, pathway, target, and bioactivity information.

**Key columns in `final_candidates_info.csv`:**

| Column | Description |
|--------|-------------|
| `seed_name` | Name of the matched seed compound |
| `HepG2-FFA模型证据` | Literature evidence for lipid-lowering in HepG2-FFA model |
| `targetmol_id` | TargetMol compound ID for submission |
| `similarity` | ECFP4 Tanimoto similarity to seed |
| `hepg2_toxic_proba_S2DV` | HepG2 toxicity probability (S2DV model) |
| `hepg2_toxic_proba_nature` | HepG2 toxicity probability (Nature chemprop model) |
| `Name`, `Pathways`, `Target`, `Bioactivity` | Agent-curated annotations for candidate compounds |

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

| Model | Repository | Reference |
|-------|------------|-------|
| S2DV (HepG2 cytotoxicity) | [NTU-MedAI/S2DV](https://github.com/NTU-MedAI/S2DV) | Shao et al., *Briefings in Bioinformatics*, 2022 |
| Chemprop cytotoxicity | [felixjwong/antibioticsai](https://github.com/felixjwong/antibioticsai) | Wong et al., *Nature*, 2024 |
| Compound library | TargetMol T001 (~22,966 compounds) | — |

---

## License

This repository is for competition submission and academic reference only.
The pre-trained model weights in `S2DV-main/` and `antibioticsai/` are subject to their respective upstream licenses.
