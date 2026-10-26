# -*- coding: utf-8 -*-
import os
import re

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SUPP = os.path.join(BASE, "paper", "v51_ReliFuse_Supplementary.md")

with open(SUPP, "r", encoding="utf-8") as f:
    p = f.read()

H = "## "

def note(heading_number, text):
    # heading_number is ASCII digit; the full heading starts with "## 补充表 N."
    prefix = H + "\u8865\u5145\u8868 " + heading_number + r"\."
    pat = "(?m)^(" + prefix + ".*)$"
    global p
    p = re.sub(pat, lambda m: m.group(1) + text, p)

note("3", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 IX\uff1b\u6b64\u5904\u4e3a\u5b8c\u6574\u660e\u7ec6\uff09")
note("9", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 XIV\uff09")
note("10", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 XIV\uff09")
note("13", "\uff08\u6b63\u6587\u8868 XIII \u4e3a\u7cbe\u70bc\u6458\u8981\uff09")
note("16", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 III\uff1b\u6b64\u5904\u4e3a\u5b8c\u6574\u660e\u7ec6\uff09")
note("20", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 X\uff1b\u6b64\u5904\u4e3a\u5b8c\u6574\u660e\u7ec6\uff09")
note("24", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 IV\uff1b\u6b64\u5904\u4e3a\u5b8c\u6574\u660e\u7ec6\uff09")

with open(SUPP, "w", encoding="utf-8") as f:
    f.write(p)

print("done")
