# Code

This folder regenerates the v51 ReliFuse manuscript assets from `../data/` and
`../results/`. Run scripts from this directory.

## Manuscript outputs

| Script | Output |
|---|---|
| `build_v51_docx.py` | `paper/v51_ReliFuse_paper.docx` and `paper/v51_ReliFuse_Supplementary.docx` |
| `build_latex.py` | `paper/JBHI_template/IEEE-TJ-color-latex-template/ReliFuse.tex` |
| `check_refs.py` | Main/supplementary figure, table, and cross-reference audit |
| `check_manuscript.py` | Current-v51 structural and forbidden-token audit |

The LaTeX PDF is compiled in the template directory with `pdflatex ReliFuse`
(run twice for stable references).

## Core v51 experiment and figures

| Script | Output |
|---|---|
| `make_two_axis_complementarity.py` | Fig. 3, `results/e20_two_axis_complementarity.csv` |
| `make_key_cell_robustness.py` | `results/e25_key_cell_robustness.csv` (key-cell bootstrap / permutation / LOPO) |
| `make_external_recalibration.py` | `results/e26_external_recalibration*.csv` (external cohort recalibration) |
| `make_pu_fusion_control.py` | `results/e27_pu_fusion_control*.csv` (PU fusion-head controls) |
| `make_explicit_gate_valid.py` | Explicit-gate validation analysis and `results/e24_explicit_gate_valid_summary.csv` |
| `make_reliability_ablation.py` | Reliability-conditioned gate ablations |
| `make_fig1_architecture.py` | Fig. 1 |
| `make_fig2_protocol.py` | Fig. 2 |
| `make_fig_data.py` | Main benchmark figures and legacy supporting figures |
| `make_headroom_control.py` | `results/e13_quartile_headroom.csv`, `results/e14_q4_margin_strata.csv` |
| `make_stage1_analyses.py` | Homology, calibration, and budget analyses |
| `make_decision_policy.py` | Verification-budget policy outputs |
| `make_supp_figures.py` | Supplementary figures |

Older v44–v50 figure scripts are retained as provenance; displayed v51 figure
numbers are defined by `../figures/README.md`, not by script names.

## Reproducibility convention

The public benchmark package has not yet been deposited. Its four logical
components are currently represented by:

- fixed and homology-aware splits: `../data/splits/`
- per-protein predictions and features: `../data/`
- unified metrics: `metrics.py`
- command-line evaluation: `evaluate.py`

Do not create or cite a local `../benchmark/` directory until the accepted
GitHub/Zenodo release is assembled. The manuscript intentionally says the
resource accompanies the paper release rather than claiming it is already
present in this working tree.
