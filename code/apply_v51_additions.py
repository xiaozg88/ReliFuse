# -*- coding: utf-8 -*-
import os
import re

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MAIN = os.path.join(BASE, "paper", "v51_ReliFuse_paper.md")
SUPP = os.path.join(BASE, "paper", "v51_ReliFuse_Supplementary.md")


def rd(p):
    with open(p, "r", encoding="utf-8") as f:
        return f.read()


def wr(p, s):
    with open(p, "w", encoding="utf-8") as f:
        f.write(s)


D = "\u0394"
M = "\u2212"


def cluster_para():
    return (
        "To rule out homology memorization at the train/test boundary, a homology-aware cluster "
        "cross-validation (greedy k-mer clustering over all sequences, threshold 50%, five "
        "cluster-disjoint folds, five seeds) yields 25 paired fold observations (Table @CLUSTER@). "
        "ReliFuse improves MCC over MKMP by +0.0268 (Wilcoxon p = 3.0e-07, 24/1/0) and over the "
        "simple weighted fusion by " + M + "0.0032 (p = 0.141), while the weighted fusion improves "
        "over MKMP by +0.0300 (p = 4.2e-07). Across folds the mean MCC is 0.7350 (ESM-MLP), "
        "0.7485 (MKMP), 0.7785 (weighted fusion), and 0.7753 (ReliFuse)."
    )


def cluster_table():
    return (
        "Table @CLUSTER@. Homology-aware cluster cross-validation (5 seeds \u00d7 5 folds = 25 "
        "paired observations; Wilcoxon test).\n\n"
        "| Comparison | Metric | Mean " + D + " | Win/Loss/Tie | p |\n"
        "| --- | --- | ---: | ---: | ---: |\n"
        "| ReliFuse vs MKMP | MCC | +0.0268 | 24/1/0 | 3.0e-07 |\n"
        "| ReliFuse vs MKMP | F1 | +0.0196 | 24/1/0 | 1.8e-07 |\n"
        "| ReliFuse vs MKMP | AUC | +0.0077 | 24/1/0 | 1.2e-07 |\n"
        "| ReliFuse vs MKMP | Brier | " + M + "0.0090 | 0/25/0 | 6.0e-08 |\n"
        "| ReliFuse vs MKMP | ECE | " + M + "0.0060 | 6/19/0 | 1.0e-03 |\n"
        "| ReliFuse vs weighted fusion | MCC | " + M + "0.0032 | 9/16/0 | 0.141 |\n"
        "| Weighted fusion vs MKMP | MCC | +0.0300 | 24/1/0 | 4.2e-07 |\n"
    )


def homfree_table():
    return (
        "Table @HOMFREE@. Strict homology stratification by maximum test-to-train sequence identity "
        "(20-seed pooled probabilities).\n\n"
        "| Subset | n | MKMP MCC | ESM-MLP MCC | SecProGNN MCC | Weighted MCC | ReliFuse MCC | Best single pathway | "
        + D + " vs best single [95% CI] |\n"
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |\n"
        "| All test set | 1613 | 0.7814 | 0.7819 | 0.7835 | 0.8125 | 0.8176 | SecProGNN | +0.0341 [+0.0088, +0.0581] |\n"
        "| Identity < 0.20 (homology-free) | 843 | 0.7341 | 0.6928 | 0.7266 | 0.7559 | 0.7658 | MKMP | +0.0317 ["
        + M + "0.0068, +0.0753] |\n"
        "| Identity > 0.25 (homology-rich) | 672 | 0.8314 | 0.8834 | 0.8493 | 0.8755 | 0.8723 | ESM-MLP | "
        + M + "0.0111 [" + M + "0.0447, +0.0249] |\n"
    )


def strict_table():
    return (
        "Table @STRICT@. Strict-label sensitivity (20-seed pooled probabilities).\n\n"
        "| Label subset | MKMP MCC | ReliFuse MCC | " + D + "MCC |\n"
        "| --- | ---: | ---: | ---: |\n"
        "| All proteins | 0.7814 | 0.8176 | +0.036 |\n"
        "| Excluding plasma-only | 0.7876 | 0.8289 | +0.041 |\n"
        "| Excluding intracellular | 0.7809 | 0.8165 | +0.036 |\n"
        "| Extracellular-annotated only | 0.7681 | 0.8187 | +0.051 |\n"
        "| Excluding contaminant proteins | 0.7820 | 0.8167 | +0.035 |\n"
        "| Strict secreted subset | 0.7747 | 0.8277 | +0.053 |\n"
    )


def replace_table():
    return (
        "Table @REPLACE@. Structure-pathway replacement control (same fusion protocol and ESM sequence "
        "pathway; 20 shared seeds).\n\n"
        "| Structure pathway | Single-pathway MCC | Fused MCC | " + D + "(fused " + M + " single) |\n"
        "| --- | ---: | ---: | ---: |\n"
        "| SecProGNN (baseline) | 0.7516 | 0.7811 | +0.0295 |\n"
        "| mgate | 0.7568 | 0.7831 | +0.0263 |\n"
        "| MKMP (this) | 0.7609 | 0.7838 | +0.0229 |\n"
    )


def source_table():
    return (
        "Table @SOURCE@. Structure-source consistency and cross-source gating (n = 86 balanced subset).\n\n"
        "| Metric | Value |\n"
        "| --- | ---: |\n"
        "| AlphaFold MKMP AUC | 0.9145 |\n"
        "| ESMFold MKMP AUC | 0.4938 |\n"
        "| Pearson r (logit) | 0.0631 |\n"
        "| Spearman \u03c1 (logit) | " + M + "0.2124 |\n"
        "| Sign agreement | 0.5698 |\n"
        "| Sequence pathway AUC | 0.9168 \u00b1 0.0061 |\n"
        "| Consistency-gated AUC | 0.8662 \u00b1 0.0029 |\n"
    )


def go_table():
    return (
        "Table @GO@. GO enrichment of the Top-100 high-confidence undetected proteins (Fisher exact, "
        "FDR < 0.05; top terms).\n\n"
        "| GO term | Description | fg | bg | FDR |\n"
        "| --- | --- | ---: | ---: | ---: |\n"
        "| GO:0004896 | cytokine receptor activity | 4 | 50 | 6.66e-03 |\n"
        "| GO:0008543 | fibroblast growth factor receptor signaling pathway | 4 | 43 | 7.33e-03 |\n"
        "| GO:0043410 | positive regulation of MAPK cascade | 6 | 178 | 8.06e-03 |\n"
        "| GO:0005111 | type 2 fibroblast growth factor receptor binding | 2 | 5 | 8.12e-03 |\n"
        "| GO:0008083 | growth factor activity | 5 | 139 | 1.43e-02 |\n"
        "| GO:0019221 | cytokine-mediated signaling pathway | 5 | 136 | 1.58e-02 |\n"
        "| GO:0060333 | type II interferon-mediated signaling pathway | 2 | 13 | 2.27e-02 |\n"
        "| GO:0035666 | TRIF-dependent toll-like receptor signaling pathway | 2 | 11 | 1.92e-02 |\n"
    )


def main():
    s = rd(MAIN)

    # 1) cluster CV paragraph + table in IV-C, before existing homology-free paragraph
    anchor = "On the homology-free subset"
    if anchor not in s:
        raise RuntimeError("missing IV-C anchor")
    s = s.replace(anchor, cluster_para() + "\n\n" + cluster_table() + "\n\n" + anchor, 1)

    # 2) reference the new homology-free table inside existing paragraph
    s = s.replace(
        "and ReliFuse reaches 0.7658, significantly higher",
        "and ReliFuse reaches 0.7658 (Table @HOMFREE@), significantly higher",
        1,
    )
    # place the homology-free table after the IV-C paragraph
    s = s.replace(
        "does not come from homology memorization.",
        "does not come from homology memorization.\n\n" + homfree_table(),
        1,
    )

    # 3) IV-F strict-label + replacement + inference paragraph
    start = "A strict-label sensitivity analysis"
    end = "The fusion module requires only approximately"
    si = s.find(start)
    ei = s.find(end)
    if si == -1 or ei == -1 or ei <= si:
        raise RuntimeError("missing IV-F block anchors")
    new_ivf = (
        "A strict-label sensitivity analysis (removing plasma-only, intracellular, and contaminant "
        "proteins) shows ReliFuse's advantage over MKMP is stably retained (+0.036~+0.053) across "
        "label subsets (Table @STRICT@).\n\n"
        + strict_table()
        + "\n\nReplacing MKMP with the SecProGNN baseline under the same fusion protocol and ESM "
        "pathway (Table @REPLACE@) gives SecProGNN+ESM MCC 0.7811 versus MKMP+ESM 0.7838 (paired "
        + D + " = +0.0028, p = 0.43); the mgate variant gives 0.7831 (paired " + D + " = +0.0007, "
        "p = 0.70), confirming that the framework's performance comes mainly from fusion and the "
        "sequence pathway, insensitive to the structure-pathway instance.\n\n"
        + replace_table()
        + "\n\n"
    )
    s = s[:si] + new_ivf + s[ei:]

    # 4) GO enrichment paragraph + table in Discussion A
    go_anchor = "### B. Scope and Limitations"
    if go_anchor not in s:
        raise RuntimeError("missing Discussion B anchor")
    go_para = (
        "The Top-100 high-confidence undetected proteins are enriched for secretion- and "
        "signaling-related GO terms, including cytokine receptor activity (GO:0004896), fibroblast "
        "growth factor receptor signaling (GO:0008543), and growth factor activity (GO:0008083) "
        "(Table @GO@). This indicates that the ranked \"undetected\" candidates are biologically "
        "secretion/signaling-biased leads rather than random artifacts; they remain hypotheses for "
        "targeted verification, not confirmed missed detections."
    )
    s = s.replace(go_anchor, go_para + "\n\n" + go_table() + "\n\n" + go_anchor, 1)

    # 5) structure-source consistency paragraph + table in Discussion B, before PU table
    source_anchor = "Table IX. PU label-noise sensitivity"
    if source_anchor not in s:
        raise RuntimeError("missing PU table anchor")
    source_para = (
        "Because the structure pathway depends on the structure-prediction source and feature "
        "protocol, we add a structure-source consistency gate (Section III-E). On an 86-protein "
        "balanced subset with both AlphaFold and ESMFold structures, MKMP AUC drops from 0.9145 "
        "(AlphaFold) to 0.4938 (ESMFold), with low logit correlation (Pearson r = 0.0631; Spearman "
        "\u03c1 = " + M + "0.2124) and sign agreement 0.5698; a consistency gate that uses the structure "
        "logit only when both sources agree in direction, otherwise falling back to the sequence "
        "logit, recovers AUC to 0.8662 \u00b1 0.0029, close to the sequence-only 0.9168 \u00b1 0.0061 "
        "(Table @SOURCE@)."
    )
    s = s.replace(source_anchor, source_para + "\n\n" + source_table() + "\n\n" + source_anchor, 1)

    # 6) renumber existing Roman table labels/captions/references
    mapping = {
        "IX": "XV",
        "VIII": "XII",
        "VII": "XI",
        "VI": "VIII",
        "V": "VII",
        "IV": "VI",
        "III": "V",
        "II": "II",
        "I": "I",
    }
    order = ["IX", "VIII", "VII", "VI", "IV", "III", "II", "I", "V"]
    for old in order:
        new = mapping[old]
        if new == old:
            continue
        pat = "Table " + old + "(?![A-Za-z])"
        s = re.sub(pat, "Table " + new, s)

    # 7) assign final numbers to new tables
    tokens = {
        "@CLUSTER@": "III",
        "@HOMFREE@": "IV",
        "@STRICT@": "IX",
        "@REPLACE@": "X",
        "@GO@": "XIII",
        "@SOURCE@": "XIV",
    }
    for tok, num in tokens.items():
        s = s.replace(tok, num)

    wr(MAIN, s)

    # 8) update supplementary annotations
    p = rd(SUPP)

    def rep(old, new):
        nonlocal p
        p = p.replace(old, new)

    # stale promoted annotations due to renumbering
    rep("\u8865\u5145\u8868 5. \u95e8\u63a7\u53ef\u8bc6\u522b\u6027\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 IV\uff0c\u56fe\u793a\u89c1\u8865\u5145\u56fe 8\uff09",
        "\u8865\u5145\u8868 5. \u95e8\u63a7\u53ef\u8bc6\u522b\u6027\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 VI\uff0c\u56fe\u793a\u89c1\u8865\u5145\u56fe 8\uff09")
    rep("\uff08\u5df2\u5408\u5e76\u81f3\u6b63\u6587\u8868 III\uff09",
        "\uff08\u5df2\u5408\u5e76\u81f3\u6b63\u6587\u8868 V\uff09")
    rep("\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 IX\uff1b\u8bc4\u4f30\u5c42\u9762",
        "\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 XV\uff1b\u8bc4\u4f30\u5c42\u9762")
    rep("\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 IV\u2013V",
        "\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 VI\u2013VII")

    # new promotions: append notes to headings
    def append_note(heading_prefix, note):
        nonlocal p
        pat = "(?m)^(" + heading_prefix + ".*)$"
        p = re.sub(pat, lambda m: m.group(1) + note, p)

    append_note("\u8865\u5145\u8868 3\\.", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 IX\uff1b\u6b64\u5904\u4e3a\u5b8c\u6574\u660e\u7ec6\uff09")
    append_note("\u8865\u5145\u8868 9\\.", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 XIV\uff09")
    append_note("\u8865\u5145\u8868 10\\.", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 XIV\uff09")
    append_note("\u8865\u5145\u8868 13\\.", "\uff08\u6b63\u6587\u8868 XIII \u4e3a\u7cbe\u70bc\u6458\u8981\uff09")
    append_note("\u8865\u5145\u8868 16\\.", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 III\uff1b\u6b64\u5904\u4e3a\u5b8c\u6574\u660e\u7ec6\uff09")
    append_note("\u8865\u5145\u8868 20\\.", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 X\uff1b\u6b64\u5904\u4e3a\u5b8c\u6574\u660e\u7ec6\uff09")
    append_note("\u8865\u5145\u8868 24\\.", "\uff08\u5df2\u63d0\u5347\u81f3\u6b63\u6587\u8868 IV\uff1b\u6b64\u5904\u4e3a\u5b8c\u6574\u660e\u7ec6\uff09")

    wr(SUPP, p)

    print("main tables (captions in first-appearance order):")
    for m in re.finditer(r"^Table ([IVX]+)\.", s, flags=re.M):
        print("  ", m.group(1))
    print("new tokens left:", re.findall(r"@[A-Z]+@", s))
    print("supp heading count:", len(re.findall(r"^## \u8865\u5145\u8868", p, flags=re.M)))


if __name__ == "__main__":
    main()
