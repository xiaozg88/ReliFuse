# -*- coding: utf-8 -*-
import os

# One-off patch (already applied): removes the margin element from the OMML
# gate-input equation. Kept for provenance; operates on the vendored copy.
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "omml_helpers.py")
with open(path, "r", encoding="utf-8") as f:
    s = f.read()

needle = "r('margin'), "
if needle not in s:
    raise RuntimeError("margin element not found")
s = s.replace(needle, "")

with open(path, "w", encoding="utf-8") as f:
    f.write(s)

print("done")
