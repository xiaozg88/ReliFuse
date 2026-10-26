# Data dictionary

All files are protein-level tables for the 8,069-protein BALF benchmark or for
the supporting analyses. Unless noted, rows are keyed by UniProt accession in
the `name` column and the test split contains 1,613 proteins (482 positives).

## 1. Primary per-protein predictions

| File(s) | Rows | Columns | Meaning |
|---|---|---|---|
| `preds_esm150m_v4_s<seed>.csv` (20 files, seeds 3407…1000) | 1,613 | `name, y, z_mkmp_cal, z_esm_cal, p_mkmp, p_esm, p_weighted, p_stack, p_gated, p_gated2, p_stack_iso, g_gated2, b_gated2` | main-pipeline predictions for one random seed: temperature-calibrated logits, probabilities of the structure path (`p_mkmp`), sequence path (`p_esm`), weighted fusion, logistic stacking (`p_stack`), isotonic stacking (`p_stack_iso`), and ReliFuse (`p_gated2`) together with its per-protein gate `g` and bias `b` |
| `base_mgate_preds_s<seed>.csv` (20 files) | 1,613 | `name, y, p_base, p_mgate` | SecProGNN baseline (`p_base`) and the `mgate` variant (`p_mgate`) on the same test split |
| `benchmark_pred_<method>.csv` (7 files) | 1,613 | `name, y, p` | 20-seed probability-aggregated predictions used by the leaderboard (`SecProGNN, MKMP, ESM_MLP, weighted, logistic_stacking, isotonic_stacking, ReliFuse`) |
| `esm150m_v4_fusion_results.csv` | 140 | per-seed metrics | per-seed discrimination/calibration/utility metrics for 7 methods (main tables 2, 5, 6) |
| `clustercv_esm150m_fusion_results.csv` | 20 | per-fold metrics | homology-aware cluster CV, 5 folds × 4 methods (main table 9) |

## 2. Structural features and mechanism analyses

| File(s) | Rows | Meaning |
|---|---|---|
| `plddt_features_8069.npy`, `plddt_names_8069.npy` | 8,069 | per-protein AlphaFold pLDDT summary (mean, fraction low/medium/high) and the matching accession list |
| `trajectory_v12_s<seed>.csv` (5 files) | 1,613 | per-protein, per-layer, per-kernel gate weights `a<layer>_<kernel>` (3 layers × 4 kernels) of the MKMP structure path |
| `mkmp_v12_fingerprint.csv` | 1,613 | biophysical fingerprint (pLDDT, length, hydrophobicity, aromaticity, entropy, N-terminal hydrophobicity, distance moments) plus the mechanism-divergence score `D` |
| `mkmp_v12_bio_corr.csv` | 24 | correlation between each kernel weight and each biophysical feature (Figure 6c) |
| `homology_identity.csv` | 1,613 | maximum sequence identity of each test protein to the training set (MMseqs2) |

## 3. Supporting analyses

| File | Meaning |
|---|---|
| `e1_fair_transfer_results.csv` | GlycoHBF fairness-transfer protocols (zero-shot / fine-tuned / retrained) |
| `e2_overlap_strata_metrics.csv` | GlycoHBF overlap/homology strata (supplementary Table 4) |
| `e8_label_sensitivity_metrics.csv` | strict-label sensitivity strata (supplementary Table 3) |
| `e9_shap_mean_abs.csv` | mean |SHAP| of the 16 gate features |
| `e10_counterfactual.csv` | counterfactual interventions on the gate (pLDDT injection, path swap, disagreement removal) |
| `e12_gate_identifiability.csv` | fixed/permuted gate ablations (supplementary Figure 8) |
| `e13_gate_path_reliability_quintiles.csv` | gate quintiles vs sequence-path advantage |
| `structure_safety_switch.csv` | AlphaFold vs ESMFold structure-source consistency (5 seeds) |
| `dca_netbenefit.csv` | decision-curve net benefit by threshold |
| `peptide_pilot_103.csv` | per-peptide Prosit iRT / fragment intensity / spectral quality / proteotypicity for the 103-pilot set |

## 4. Not shipped here

| Artifact | Size | Note |
|---|---|---|
| `data_list.pt` (AlphaFold contact graphs, 8,069 proteins) | ≈1.05 GB | see `DOWNLOAD.md` |
| `esm_meanpool_8069_esm2_150m.npz`, `..._esm2_650m_faesm.npz`, `..._esmc300m.npz` | ≈90 MB | see `DOWNLOAD.md`; fully regenerable from UniProt sequences |

## 5. Provenance

Detection labels and protein sets derive from HBFP and the SecProGNN study;
structures from AlphaFold DB; the cross-fluid companion benchmark (GlycoHBF) and
the two external BALF cohorts (PMC5357681, PXD012645) are described in the manuscript Data and Code Availability section; the public
data card will accompany the accepted release. Please cite the primary resources
when redistributing.
