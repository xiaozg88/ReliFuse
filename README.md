# ReliFuse — Reproducibility Resource

Benchmark resource accompanying the manuscript:

> **ReliFuse: Reliability-Gated Structure–Sequence Complementarity for
> Budget-Constrained Protein Candidate Ranking in BALF Mass Spectrometry**
> (submitted to *IEEE Journal of Biomedical and Health Informatics*)

ReliFuse frames BALF mass-spectrometry protein detection as a positive-unlabeled
(PU), budget-constrained candidate-ranking problem, and shows that structure and
sequence pathways are complementary along two orthogonal reliability axes —
structural confidence (pLDDT) and evolutionary support (sequence identity) —
with the fusion gain concentrated in the low-pLDDT / low-identity cell and
realized by a minimal calibrated operator (temperature scaling + per-protein
bias; approximately 3.2K trainable parameters).

## Repository layout

| Folder | Content |
|---|---|
| `code/` | All table/figure-generation, evaluation, and QA scripts (see `code/README.md`) |
| `data/` | Fixed and homology-aware splits, 20-seed per-protein predictions, pLDDT and homology features, external-cohort predictions |
| `results/` | All experiment CSVs (E1–E27), leaderboard, logits caches |
| `figures/` | All manuscript and supplementary figures (PNG/PDF) |
| `DATA_CARD.md` | Dataset provenance, label semantics (PU disclosure), limitations |
| `LICENSE` | MIT |
| `requirements.txt` | Python dependencies |

## Quick start

```bash
pip install -r requirements.txt
cd code

# Core v51 analyses
python make_two_axis_complementarity.py   # Fig. 3 + results/e20_*.csv
python make_key_cell_robustness.py        # Supp. Table 32 (~1 min)
python make_external_recalibration.py     # Supp. Table 33 (~25 s)
python make_pu_fusion_control.py          # Supp. Table 34 (~5 min)

# Manuscript QA
python check_manuscript.py                # FAIL=0 WARN=0 OK=24 expected
python check_supplement_en.py             # FAIL=0 expected
```

## Reporting conventions (important)

Two metric conventions coexist and must not be mixed:

- **Seed-wise means** over 20 shared random seeds — used in the main tables
  (e.g., MCC 0.8019 / AUC 0.9687 / ECE 0.0184);
- **20-seed pooled probabilities** — used for ROC/PR curves, the leaderboard,
  Top-K budget analysis, and homology stratification
  (e.g., MCC 0.8176 / AUC 0.9724 / ECE 0.0158).

All baselines reported in the manuscript's Table I (DeepSec, ESMSec, GraphIDN,
ESM-2 650M + MLP, ESM-C 300M + MLP, GAT, GIN, SecProGNN) were trained in-house
on this benchmark under the same split and metric protocol.

## Scope and honest boundaries

ReliFuse is a candidate-prioritization tool **within the same BALF workflow**.
Direct cross-cohort/cross-fluid deployment is not supported without
recalibration (see DATA_CARD.md and the manuscript's Section V).

## Data availability

Large/raw artifacts (upstream pipeline outputs, per-seed model weights) are
archived on Zenodo: `<DOI to be inserted on acceptance>`. See
`data/DOWNLOAD.md`.

## License

MIT (see `LICENSE`).
