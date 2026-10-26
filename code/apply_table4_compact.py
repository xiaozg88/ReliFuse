# -*- coding: utf-8 -*-
import os

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MAIN = os.path.join(BASE, "paper", "v51_ReliFuse_paper.md")

with open(MAIN, "r", encoding="utf-8") as f:
    s = f.read()

old_header = (
    "| Subset | n | MKMP MCC | ESM-MLP MCC | SecProGNN MCC | Weighted MCC | "
    "ReliFuse MCC | Best single pathway | \u0394 vs best single [95% CI] |"
)
new_header = (
    "| Subset | n | MKMP | ESM-MLP | SecProGNN | Weighted | ReliFuse | "
    "Best single | \u0394 vs best [95% CI] |"
)
if old_header not in s:
    raise RuntimeError("old table IV header not found")
s = s.replace(old_header, new_header, 1)

s = s.replace("| All test set | 1613 |", "| All | 1613 |", 1)
s = s.replace("| Identity < 0.20 (homology-free) | 843 |", "| < 0.20 (homology-free) | 843 |", 1)
s = s.replace("| Identity > 0.25 (homology-rich) | 672 |", "| > 0.25 (homology-rich) | 672 |", 1)

old_cap = "Table IV. Strict homology stratification by maximum test-to-train sequence identity (20-seed pooled probabilities)."
new_cap = "Table IV. Strict homology stratification by maximum test-to-train sequence identity (20-seed pooled MCC)."
if old_cap not in s:
    raise RuntimeError("old table IV caption not found")
s = s.replace(old_cap, new_cap, 1)

with open(MAIN, "w", encoding="utf-8") as f:
    f.write(s)

print("done")
