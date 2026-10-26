# -*- coding: utf-8 -*-
import os

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MAIN = os.path.join(BASE, "paper", "v51_ReliFuse_paper.md")

with open(MAIN, "r", encoding="utf-8") as f:
    s = f.read()

old_para = (
    "Cross-cohort and cross-fluid comparisons delimit the tool's scope rather than claim "
    "generalization. On two independent BALF cohorts never involved in training, the sequence "
    "pathway alone achieves AUC 0.845 (Hmmier 2017 [34]) and 0.668 (PXD012645 [33]), whereas "
    "the full fusion model reaches only 0.520 and 0.591\u2014below the sequence-only pathway\u2014"
    "because the structure pathway is unstable under cross-platform structure-feature protocols "
    "(Table XII). On GlycoHBF (15 non-BALF fluids), direct transfer yields AUC \u2248 0.50. "
    "ReliFuse is therefore a candidate-prioritization tool within the same BALF pipeline; its "
    "parameters must be refit per pipeline, and direct cross-cohort deployment is not supported."
)

new_para = (
    "Cross-cohort and cross-fluid comparisons delimit the tool's scope rather than claim "
    "generalization. On two independent BALF cohorts never involved in training, the sequence "
    "pathway alone achieves AUC 0.845 (Hmmier 2017 [34]) and 0.668 (PXD012645 [33]), whereas "
    "the full fusion model reaches only 0.520 and 0.591\u2014below the sequence-only pathway "
    "(Table XII). On GlycoHBF (15 non-BALF fluids), direct transfer yields AUC \u2248 0.50 "
    "(Table XII). The structure-pathway fragility that drives these cross-source and "
    "cross-cohort drops is quantified in Table XIV. ReliFuse is therefore a candidate-"
    "prioritization tool within the same BALF pipeline; its parameters must be refit per "
    "pipeline, and direct cross-cohort deployment is not supported."
)

if old_para not in s:
    raise RuntimeError("old cross-cohort paragraph not found")
s = s.replace(old_para, new_para, 1)

old_caption = "Table XII. Structure-source and cross-cohort/cross-fluid out-of-domain AUC."
new_caption = "Table XII. Cross-cohort/cross-fluid out-of-domain AUC."
if old_caption not in s:
    raise RuntimeError("old Table XII caption not found")
s = s.replace(old_caption, new_caption, 1)

old_row = "| ESMFold (cross-source) | 0.4938 | \u2013 | \u2013 | 0.8662 |\n"
if old_row not in s:
    raise RuntimeError("old ESMFold row not found")
s = s.replace(old_row, "", 1)

with open(MAIN, "w", encoding="utf-8") as f:
    f.write(s)

print("done")
