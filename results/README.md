# Result tables

Experiment outputs backing the v51 manuscript. Files are CSV unless otherwise
noted. Table numbers refer to the current v51 manuscript.

## v51 core analyses

| File | Manuscript object |
|---|---|
| `e20_two_axis_complementarity.csv` | Fig. 3 / Supplementary Table 0 — pLDDT × identity complementarity map |
| `e20_axis_orthogonality.csv` | Section IV-A — Spearman/Pearson correlation between the two reliability axes (pLDDT × max identity) on the test set |
| `e25_key_cell_robustness.csv` | Supplementary Table 32 — key-cell bootstrap CI / permutation / LOPO robustness |
| `e26_external_recalibration.csv`, `e26_external_recalibration_paired.csv` | Supplementary Table 33 — external-BALF pathway recalibration (feature failure vs scale drift) |
| `e27_pu_fusion_control.csv`, `e27_pu_fusion_control_summary.csv` | Supplementary Table 34 — PU fusion-head controls (uPU/nnPU/bagging/two-step) |
| `e21_two_by_two_robustness.csv` | Supplementary Table 30 — coarse 2×2 robustness check |
| `e22_gate_correlations.csv` | Supplementary Table 31 — gate correlations with pLDDT and identity |
| `e23_explicit_gate_summary.csv` | Explicit-gate ablation |
| `e24_explicit_gate_valid_summary.csv` | Table VII / Supplementary Table 31 — validation-fitted explicit gate |
| `val_logits_cache.csv` | Per-protein validation logits used to fit explicit gates |
| `leaderboard.csv` | Supplementary Table 26 — 20-seed aggregate leaderboard |

## Main benchmark and controls

| File | Manuscript object |
|---|---|
| `e1_plddt_stratification.csv` / `.md` | Supplementary Table 14 — pLDDT-stratified metrics |
| `e2_paired_stats.csv` | Supplementary Table 15 — paired Wilcoxon/Holm statistics |
| `e3_clustercv_stats.csv` | Homology-aware cluster-CV statistics |
| `e4_localization_shortcut.csv` | Supplementary Table 18 — localization-shortcut control |
| `e5_topk_budget.csv`, `e5b_budget_gain.csv` | Budget analyses |
| `e6_pu_sensitivity.csv` | Table XV / Supplementary Table 19 — PU label-noise sensitivity |
| `e7_substitution_control*.csv` / `.md` | Table X / Supplementary Table 20 — structure-path replacement |
| `e8_multiseed_clustercv*.csv` | Table III / Supplementary Table 16 — 5 seed × 5 fold cluster CV |
| `e9_pxd012645_full_model.csv` | Supplementary Table 8 — external BALF cohort |
| `e13_quartile_headroom.csv` | Supplementary Table 21 — disagreement/headroom control |
| `e14_q4_margin_strata.csv` | Supplementary Table 22 — decision-boundary margin control |
| `e15_budget_utility.csv` | Table V and budget-utility analyses |
| `e16_homology_free.csv` | Table IV / Supplementary Table 24 — strict homology strata |
| `e17_calibration_headtohead.csv`, `e17_calibration_paired.csv` | Supplementary Table 25 — fusion head-to-head |
| `e18_planning_accuracy.csv` | Supplementary Table 27 — yield-planning error |
| `e19_budget_optimality.csv` | Supplementary Table 28 — optimal verification depth |

`e7_secpognn_control.csv`, `e7_summary.csv`, `v45_experiment_report.md`, and
`supp_text.txt` are historical QA/provenance records and should not be treated
as current v51 manuscript tables.

Regenerate data-derived outputs with the scripts in `../code/` (see
`../code/README.md`).
