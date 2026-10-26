# -*- coding: utf-8 -*-
import os

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MAIN = os.path.join(BASE, "paper", "v51_ReliFuse_paper.md")

with open(MAIN, "r", encoding="utf-8") as f:
    s = f.read()

reps = [
    ("| Metric | Learned gate | Explicit gate (pLDDT + identity) |",
     "| Metric | Learned gate | Explicit gate |"),
    ("| Key-cell \u0394MCC (<70 \u00d7 <0.20) | +0.200 | +0.114 |",
     "| Key-cell \u0394MCC | +0.200 | +0.114 |"),
    ("| Extracellular-annotated only | 0.7681 | 0.8187 | +0.051 |",
     "| Extracellular-only | 0.7681 | 0.8187 | +0.051 |"),
    ("| Excluding contaminant proteins | 0.7820 | 0.8167 | +0.035 |",
     "| Excluding contaminants | 0.7820 | 0.8167 | +0.035 |"),
    ("| Strict secreted subset | 0.7747 | 0.8277 | +0.053 |",
     "| Strict secreted | 0.7747 | 0.8277 | +0.053 |"),
    ("| Structure pathway | Single-pathway MCC | Fused MCC | \u0394(fused \u2212 single) |",
     "| Structure pathway | Single MCC | Fused MCC | \u0394MCC |"),
    ("| SecProGNN (baseline) | 0.7516 | 0.7811 | +0.0295 |",
     "| SecProGNN | 0.7516 | 0.7811 | +0.0295 |"),
    ("| MKMP (this) | 0.7609 | 0.7838 | +0.0229 |",
     "| MKMP | 0.7609 | 0.7838 | +0.0229 |"),
]

for old, new in reps:
    if old not in s:
        raise RuntimeError("not found: " + old)
    s = s.replace(old, new, 1)

old_note = "Note: the explicit affine gate uses only the two reliability axes (pLDDT and sequence identity) and is fitted on the validation set; its overall AUC matches the learned gate and its MCC differs by 0.002, confirming that the two axes are sufficient for the fusion's performance."
new_note = "Note: the explicit affine gate uses only the two reliability axes (pLDDT and sequence identity) and is fitted on the validation set; the key cell is pLDDT < 70 and identity < 0.20. Its overall AUC matches the learned gate and its MCC differs by 0.002, confirming that the two axes are sufficient for the fusion's performance."
if old_note not in s:
    raise RuntimeError("note not found")
s = s.replace(old_note, new_note, 1)

with open(MAIN, "w", encoding="utf-8") as f:
    f.write(s)

print("done")
