# -*- coding: utf-8 -*-
import os

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SUPP = os.path.join(BASE, "paper", "v51_ReliFuse_Supplementary.md")

with open(SUPP, "r", encoding="utf-8") as f:
    p = f.read()

p = p.replace("\u8be5\u8868\u4e0e\u6b63\u6587\u8868 III", "\u8be5\u8868\u4e0e\u6b63\u6587\u8868 V")

with open(SUPP, "w", encoding="utf-8") as f:
    f.write(p)

print("done")
