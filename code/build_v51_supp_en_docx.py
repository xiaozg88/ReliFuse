# -*- coding: utf-8 -*-
"""Build the English v51 supplementary DOCX from Markdown."""
import os
from build_v51_docx import build

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
MD = os.path.join(ROOT, "paper", "v51_ReliFuse_Supplementary_EN.md")
OUT = os.path.join(ROOT, "paper", "v51_ReliFuse_Supplementary_EN.docx")

if __name__ == "__main__":
    build(MD, OUT, title_size=14)
