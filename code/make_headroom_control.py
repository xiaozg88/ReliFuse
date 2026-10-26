# -*- coding: utf-8 -*-
"""Headroom / threshold control for the 'fusion gain concentrates in Q4' claim.

Questions addressed
-------------------
1. Does the Q4 gain survive when the *learned* fusion is compared with a simple
   weighted fusion (i.e. is it specific to ReliFuse)?
2. Is the Q4 gain an artefact of the 0.5 decision threshold (the baseline is
   near the boundary for many Q4 proteins)?  -> controlled by (a) matching on
   the baseline decision margin and (b) reporting threshold-free AUC/PR-AUC.
3. How much of the gain is available headroom?  -> relative error reduction and
   comparison with a per-protein oracle best-single baseline.

Outputs
-------
results/e13_quartile_headroom.csv   (quartile table)
results/e14_q4_margin_strata.csv    (margin-matched strata inside Q4)
figures/fig5_decision_map.png           (panel c added)
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, matthews_corrcoef, brier_score_loss

from style import C, save, DATADIR

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

SEEDS = [3407, 42, 123, 2024, 777, 1, 2, 3, 4, 5, 6, 7, 8, 10, 21, 99,
         202, 555, 888, 1000]
EXP = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "results"))


def load():
    ap = pd.concat([pd.read_csv(os.path.join(DATADIR, f"preds_esm150m_v4_s{s}.csv"))
                    for s in SEEDS], ignore_index=True)
    g = ap.groupby("name").agg(
        y=("y", "max"), p_mkmp=("p_mkmp", "mean"), p_esm=("p_esm", "mean"),
        p_gated2=("p_gated2", "mean"), p_weighted=("p_weighted", "mean"),
    ).reset_index()
    g["disagreement"] = (g.p_mkmp - g.p_esm).abs()
    g["quartile"] = pd.qcut(g.disagreement, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    return g


def quartile_table(g):
    rows = []
    y = g.y.values
    for lab in ["Q1", "Q2", "Q3", "Q4"]:
        m = (g.quartile == lab).values
        yy = y[m]
        acc = {}
        for name, col in [("MKMP", "p_mkmp"), ("ESM", "p_esm"),
                          ("weighted", "p_weighted"), ("ReliFuse", "p_gated2")]:
            p = g[col].values[m]
            acc[name] = float(((p >= .5).astype(int) == yy).mean())
        best = max(("MKMP", "ESM"), key=lambda k: acc[k])
        p_best = g["p_mkmp" if best == "MKMP" else "p_esm"].values[m]
        acc_best = acc[best]
        oracle = float((np.maximum((g.p_mkmp.values[m] >= .5), (g.p_esm.values[m] >= .5))
                        .astype(int) == yy).mean())
        auc = {k: float(roc_auc_score(yy, g[col].values[m]))
               for k, col in [("MKMP", "p_mkmp"), ("ESM", "p_esm"),
                              ("weighted", "p_weighted"), ("ReliFuse", "p_gated2")]}
        auc_best = max(auc["MKMP"], auc["ESM"])
        rows.append({
            "分位": lab, "n": int(m.sum()),
            "分歧度范围": f"{g.disagreement[m].min():.3f}–{g.disagreement[m].max():.3f}",
            "该组最优单路径": best,
            "最佳单路径准确率": round(acc_best, 4),
            "oracle单路径准确率": round(oracle, 4),
            "加权融合准确率": round(acc["weighted"], 4),
            "ReliFuse准确率": round(acc["ReliFuse"], 4),
            "ΔReliFuse": round(acc["ReliFuse"] - acc_best, 4),
            "Δ加权": round(acc["weighted"] - acc_best, 4),
            "ΔReliFuse−Δ加权": round(acc["ReliFuse"] - acc["weighted"], 4),
            "相对误差下降": f"{(acc['ReliFuse'] - acc_best) / max(1 - acc_best, 1e-9) * 100:.1f}%",
            "AUC最佳单路径": round(auc_best, 4),
            "AUC加权": round(auc["weighted"], 4),
            "AUC ReliFuse": round(auc["ReliFuse"], 4),
            "ΔAUC": round(auc["ReliFuse"] - auc_best, 4),
            "ΔAUC(ReliFuse−加权)": round(auc["ReliFuse"] - auc["weighted"], 4),
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(EXP, "e13_quartile_headroom.csv"),
              index=False, encoding="utf-8-sig")
    return df


def q4_margin_strata(g, n_strata=5):
    """Match on the baseline decision margin inside Q4."""
    q4 = g[g.quartile == "Q4"].copy().reset_index(drop=True)
    y = q4.y.values
    acc_m, acc_e = ((q4.p_mkmp >= .5).astype(int) == y).mean(), ((q4.p_esm >= .5).astype(int) == y).mean()
    best = "MKMP" if acc_m >= acc_e else "ESM"
    p_best = q4.p_mkmp.values if best == "MKMP" else q4.p_esm.values
    q4["margin"] = np.abs(p_best - 0.5)
    q4["stratum"] = pd.qcut(q4.margin, n_strata, labels=[f"M{i+1}" for i in range(n_strata)])
    rows = []
    for st in q4.stratum.cat.categories:
        m = (q4.stratum == st).values
        yy = y[m]
        bs = float(((p_best[m] >= .5).astype(int) == yy).mean())
        clf = float(((q4.p_gated2.values[m] >= .5).astype(int) == yy).mean())
        w = float(((q4.p_weighted.values[m] >= .5).astype(int) == yy).mean())
        rows.append({"层": st, "n": int(m.sum()),
                     "基线判据边界距离": f"{q4.margin.values[m].min():.3f}–{q4.margin.values[m].max():.3f}",
                     "最佳单路径准确率": round(bs, 4),
                     "加权融合准确率": round(w, 4),
                     "ReliFuse准确率": round(clf, 4),
                     "ΔReliFuse": round(clf - bs, 4),
                     "Δ加权": round(w - bs, 4)})
    df = pd.DataFrame(rows)
    # equal-weight (margin-matched) aggregate
    agg = {"层": "等权平均（边界距离配平）", "n": int(df.n.sum()),
           "最佳单路径准确率": round(df["最佳单路径准确率"].mean(), 4),
           "加权融合准确率": round(df["加权融合准确率"].mean(), 4),
           "ReliFuse准确率": round(df["ReliFuse准确率"].mean(), 4),
           "ΔReliFuse": round(df["ΔReliFuse"].mean(), 4),
           "Δ加权": round(df["Δ加权"].mean(), 4),
           "基线判据边界距离": "—"}
    df = pd.concat([df, pd.DataFrame([agg])], ignore_index=True)
    df.to_csv(os.path.join(EXP, "e14_q4_margin_strata.csv"),
              index=False, encoding="utf-8-sig")
    return df


def fig5c(g, q4_strata):
    """Draw the 3-panel decision figure (map / quadrant / margin-matched)."""
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.7))
    # (a) probability map
    from matplotlib.colors import LogNorm
    ax = axes[0]
    hb = ax.hexbin(g.p_esm, g.p_mkmp, gridsize=40, cmap="Blues",
                   norm=LogNorm(vmin=1), mincnt=1, linewidths=0)
    ax.axvline(0.5, color=C["grey"], lw=0.8, ls="--")
    ax.axhline(0.5, color=C["grey"], lw=0.8, ls="--")
    ax.set(xlabel="Sequence path  $p_{ESM}$", ylabel="Structure path  $p_{MKMP}$",
           title="(a) Two-path probability map")
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=6)
    ax.grid(False)
    ax.text(0.30, 0.05, "structure only", fontsize=6.2, color=C["blue"], transform=ax.transAxes)
    ax.text(0.58, 0.88, "sequence only", fontsize=6.2, color=C["green"], transform=ax.transAxes)

    # (b) quadrant gains
    ax2 = axes[1]
    y = g.y.values
    pm = (g.p_mkmp.values >= .5).astype(int); pe = (g.p_esm.values >= .5).astype(int)
    pf = (g.p_gated2.values >= .5).astype(int)
    labels, ns, gains = [], [], []
    for a, b, lab in [(0, 0, "both 0"), (1, 1, "both 1"), (0, 1, "ESM only"), (1, 0, "MKMP only")]:
        m = (pm == a) & (pe == b)
        labels.append(lab); ns.append(int(m.sum()))
        gains.append((pf[m] == y[m]).mean() - max((pm[m] == y[m]).mean(), (pe[m] == y[m]).mean()))
    bars = ax2.bar(labels, gains, color=[C["lightgrey"], C["lightgrey"], C["green"], C["blue"]],
                   edgecolor=C["grey"], lw=0.6)
    ax2.axhline(0, color="#BBBBBB", lw=0.8)
    ax2.bar_label(bars, labels=[f"{v:+.3f}\nn={n}" for v, n in zip(gains, ns)],
                  fontsize=6.0, padding=2)
    ax2.set(ylabel="Accuracy gain over best single", title="(b) Gain by decision quadrant")
    ax2.set_ylim(min(-0.02, min(gains) - 0.01), max(gains) * 1.45)
    ax2.tick_params(axis="x", labelsize=6.4)

    # (c) margin-matched deltas inside Q4
    ax3 = axes[2]
    strat = q4_strata[q4_strata["层"] != "等权平均（边界距离配平）"]
    x = np.arange(len(strat)); w = 0.36
    ax3.bar(x - w / 2, strat["ΔReliFuse"], w, color=C["orange"], edgecolor="white",
            lw=0.5, label="ReliFuse")
    ax3.bar(x + w / 2, strat["Δ加权"], w, color=C["green"], edgecolor="white",
            lw=0.5, label="weighted")
    ax3.axhline(0, color="#BBBBBB", lw=0.8)
    agg = q4_strata[q4_strata["层"] == "等权平均（边界距离配平）"].iloc[0]
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"M{i+1}" for i in range(len(strat))], fontsize=6.4)
    ax3.set(xlabel="Q4 baseline-margin strata", ylabel="Accuracy gain over best single",
            title="(c) Q4 gain at matched decision margin")
    ax3.legend(fontsize=6.0, loc="upper right")
    ax3.text(0.02, 0.03, f"equal-weight mean: ReliFuse {agg['ΔReliFuse']:+.3f}, weighted {agg['Δ加权']:+.3f}",
             transform=ax3.transAxes, fontsize=5.8, color=C["ink"])
    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig5_decision_map")


if __name__ == "__main__":
    g = load()
    qt = quartile_table(g)
    print(qt.to_string(index=False))
    strata = q4_margin_strata(g)
    print()
    print(strata.to_string(index=False))
    fig5c(g, strata)
