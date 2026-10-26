# DATA CARD — ReliFuse v51 Benchmark Resource

## 1. Purpose
Per-protein benchmark supporting *Reliability-Gated Structure-Sequence Complementarity
for Positive-Unlabeled Candidate Ranking of BALF-Detected Proteins* (v51). It backs every
table and figure in the manuscript and supplement (Supplementary Tables 0-34).

## 2. Source and attribution
- Main benchmark: 8,069-protein BALF dataset built by SecProGNN [1] on the HBFP
  database [2]. Cite both when reusing.
- Independent external BALF cohort: Hmmier et al. 2017, BBA Clinical 7:97–104
  (PMC5357681; main-text ref. [34]);
  per-protein pathway logits in `04_data/external/bal_ext_predictions.csv`.
- Second external cohort (PXD012645, PRIDE) and GlycoHBF are referenced in the
  supplement; their per-protein predictions are not part of this package
  (aggregate results only, `05_results/e9_pxd012645_full_model.csv`).

## 3. Label semantics (IMPORTANT)
- Positive (label 1) = detected by the source BALF mass-spectrometry workflow.
- Negative (label 0) = NOT DETECTED by that workflow (detection limit, peptide
  coverage, or search-engine artefacts). These are UNLABELED, not biological
  negatives. The task is positive-unlabeled (PU) candidate ranking.
- Do not interpret model output as absolute detection probability outside the
  source workflow; recalibrate per pipeline (Supplementary Table 33).

## 4. Splits
- `04_data/splits/main_split_seed3407.csv` - fixed split, seed 3407:
  train 4,839 / validation 1,615 / test 1,613 (labels: 1,451 / 474 / 482).
  The split file lists 8,067 of the 8,069 dataset entries: two accessions
  (P01562, Q3SXR2) lack retrievable sequences and are excluded; all 2,407
  positives are retained.
- `04_data/splits/cluster_folds_50pct.csv` - homology-aware 50% k-mer clusters,
  5 cluster-disjoint folds.
- `04_data/splits/homology_identity_test.csv` - MMseqs2 max test-to-train identity.
- Test set is evaluated exactly once per model; temperatures and gates are
  validation-fitted.

## 5. Per-protein predictions (20 shared seeds, pooled)
`benchmark_pred_{SecProGNN,MKMP,ESM_MLP,weighted,logistic_stacking,isotonic_stacking,ReliFuse}.csv`
columns: `name` (UniProt accession), `y` (observational label), `p` (pooled probability).

## 6. Features
- `plddt_features_8069.npy` / `plddt_names_8069.npy` - AlphaFold pLDDT features
  (column 0 = mean pLDDT; cols 1-3 = fraction of residues with pLDDT > 50/70/90).
- `homology_identity.csv` - max sequence identity of each protein to the training set.

## 7. Aggregate experiment outputs
`05_results/e1..e27*.csv` map one-to-one onto manuscript tables; see
`05_results/README.md` (copied from the working tree) for the table-to-file mapping.
New in this revision:
- `e25_key_cell_robustness.csv` - Supplementary Table 32 (bootstrap / permutation / LOPO).
- `e26_external_recalibration*.csv` - Supplementary Table 33 (external recalibration).
- `e27_pu_fusion_control*.csv` - Supplementary Table 34 (uPU / nnPU / bagging / two-step).

## 8. Honest limitations carried by this card
- PU labels: 10-30% hidden-positive sensitivity in Supplementary Table 19 / Table XV.
- Homology: cluster CV + homology-free subset controls (Tables III-IV).
- Cross-cohort/cross-fluid negative transfer (Table XII, Supp. Tables 33).
- Structure-source fragility (ESMFold, Table XIV).
- Standard PU objectives (uPU/nnPU) do not improve the fusion head here (Supp. Table 34).

## 9. References
[1] Shao et al., SecProGNN, IEEE JBHI 29(8):6114-6123, 2025.
[2] Shao et al., HBFP, Database 2021:baab065, 2021.
