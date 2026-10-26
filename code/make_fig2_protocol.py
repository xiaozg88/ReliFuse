# -*- coding: utf-8 -*-
"""Figure 2 (v46): data, training and evaluation protocol schematic."""

import matplotlib.pyplot as plt

from style import C, rbox, arrow, save, blank_axes

STAGES = [
    ("1   Data & benchmarks", C["sky"], "#EAF4FB", [
        "Main BALF benchmark — 8,069 proteins\n(2,407 positive / 5,662 negative); split 4,841 / 1,615 / 1,613",
        "Homology-aware splits — k-mer clustering 50%,\n5-fold cluster-disjoint; MMseqs2 identity strata",
        "External & cross-fluid — GlycoHBF\n(15 body fluids, 6,938); 2 independent BALF cohorts",
    ]),
    ("2   Modeling & training", C["blue"], "#E9F2F9", [
        "Inputs — residue contact graph\n(Cα–Cα < 8.0 Å) and protein sequence",
        "Two branches — trained MKMP structure\npath + frozen PLM with trainable head",
        "Protocol — Adam 1e-3, batch 64, dropout 0.2,\nearly stop on validation; 20 shared seeds",
    ]),
    ("3   Evaluation & decision", C["green"], "#E8F6F1", [
        "Metrics — AUC / PR-AUC / MCC / F1 /\nBrier / ECE / Top-K / fixed-sensitivity / DCA",
        "Statistics — paired Wilcoxon + Holm;\n20-seed win/lose; 1,000× bootstrap",
        "Decision — Top-K budget → Prosit in-silico\npre-screen → PRM/SRM peptide list",
    ]),
]


def build():
    fig = plt.figure(figsize=(7.4, 3.9))
    ax = blank_axes(fig, xlim=(0, 1), ylim=(0, 1))

    x0s = [0.015, 0.345, 0.675]
    w = 0.31
    top, bot = 0.955, 0.335
    for (title, col, tint, cards), x0 in zip(STAGES, x0s):
        rbox(ax, x0, bot, w, top - bot, "", fc=tint, ec=col, lw=1.0)
        rbox(ax, x0, top - 0.075, w, 0.075, title, fc=col, ec=col,
             fs=8.2, bold=True, tc="white")
        ny = len(cards)
        card_h = 0.135
        gap = 0.030
        y = top - 0.075 - gap - card_h
        for text in cards:
            rbox(ax, x0 + 0.012, y, w - 0.024, card_h, text,
                 fc="white", ec=col, lw=0.7, fs=6.3)
            y -= (card_h + gap)

    # stage connectors
    arrow(ax, (0.327, 0.645), (0.343, 0.645), C["grey"], lw=1.2)
    arrow(ax, (0.657, 0.645), (0.673, 0.645), C["grey"], lw=1.2)

    # validation-only callout
    rbox(ax, 0.345, 0.225, 0.31, 0.085,
         "Validation-only fitting\ntemperature $T_1,T_2$ and the fusion gate are fit\n"
         "on validation data; the test set is used once",
         fc="white", ec=C["orange"], dashed=True, fs=6.0, tc="#7A3A00")
    arrow(ax, (0.50, 0.335), (0.50, 0.312), C["orange"], lw=0.9, ls="dashed")

    # robustness band
    rbox(ax, 0.015, 0.03, 0.97, 0.155, "", fc="#F4F4F4", ec=C["grey"], lw=0.9)
    ax.text(0.028, 0.160, "Robustness & external validation",
            fontsize=7.6, fontweight="bold", va="center", color=C["ink"])
    chips = ["Cluster-disjoint CV", "Strict label sensitivity",
             "PU label-noise ρ", "Structure-source consistency",
             "Cross-fluid & external cohorts"]
    cx = 0.028
    cw = 0.186
    for c in chips:
        rbox(ax, cx, 0.045, cw, 0.075, c, fc="white", ec=C["grey"], fs=6.2)
        cx += cw + 0.008
    arrow(ax, (0.83, 0.335), (0.83, 0.19), C["grey"], lw=1.2)

    save(fig, "fig2_protocol")


if __name__ == "__main__":
    build()
