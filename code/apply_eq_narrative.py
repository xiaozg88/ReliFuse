# -*- coding: utf-8 -*-
import os

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MAIN = os.path.join(BASE, "paper", "v51_ReliFuse_paper.md")

with open(MAIN, "r", encoding="utf-8") as f:
    s = f.read()

old_line1 = (
    "The gate input concatenates calibrated logits, probabilities, binary entropy, "
    "pathway agreement and disagreement, decision margins, structural quality, protein "
    "length, and graph density:"
)
new_line1 = (
    "The gate input concatenates calibrated logits, probabilities, binary entropy, "
    "pathway agreement and disagreement, structural quality, protein length, and graph "
    "density:"
)
if old_line1 not in s:
    raise RuntimeError("line1 not found")
s = s.replace(old_line1, new_line1, 1)

old_margin = "where margin_m = |p_m \u2212 0.5| denotes the decision margin of pathway m.\n\n"
if old_margin not in s:
    raise RuntimeError("margin line not found")
s = s.replace(old_margin, "", 1)

old_line2 = (
    "where z = (z1c, z2c), p = (p1, p2), H = (H(p1), H(p2)), and margin = (margin1, margin2) "
    "are the two pathways' calibrated logits, probabilities, binary entropies, and decision "
    "margins; s_struct = [pLDDT, f_low, f_med, f_high] packs the structural-quality features "
    "(mean pLDDT and the fractions of residues with pLDDT above 50, 70, and 90); and the "
    "protein length log(1+L) and graph density \u03c1 are also included."
)
new_line2 = (
    "where z = (z1c, z2c), p = (p1, p2), and H = (H(p1), H(p2)) are the two pathways' "
    "calibrated logits, probabilities, and binary entropies; margin = (margin1, margin2), "
    "with margin_m = |p_m \u2212 0.5|, denotes the two pathways' decision margins; s_struct = "
    "[pLDDT, f_low, f_med, f_high] packs the structural-quality features (mean pLDDT and the "
    "fractions of residues with pLDDT above 50, 70, and 90); and the protein length log(1+L) "
    "and graph density \u03c1 are also included."
)
if old_line2 not in s:
    raise RuntimeError("line2 not found")
s = s.replace(old_line2, new_line2, 1)

with open(MAIN, "w", encoding="utf-8") as f:
    f.write(s)

print("done")
