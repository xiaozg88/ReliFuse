# -*- coding: utf-8 -*-
"""Stage-2 decision policy: calibration-aware budget planning for candidate verification.

E18  planning accuracy  — model's own expected yield vs realised yield at a given budget
E19  budget optimality  — verification depth K* under a cost/value model, per method
Fig 11                  — decision-policy panels
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from style import C, save, DATADIR

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

SEEDS = [3407, 42, 123, 2024, 777, 1, 2, 3, 4, 5, 6, 7, 8, 10, 21, 99,
         202, 555, 888, 1000]
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, "..", "results"))


def load():
    ap = pd.concat([pd.read_csv(os.path.join(DATADIR, f"preds_esm150m_v4_s{s}.csv"))
                    for s in SEEDS], ignore_index=True)
    agg = ap.groupby("name").agg(
        y=("y", "max"), p_mkmp=("p_mkmp", "mean"), p_esm=("p_esm", "mean"),
        p_weighted=("p_weighted", "mean"), p_gated2=("p_gated2", "mean"),
    ).reset_index()
    bg = pd.concat([pd.read_csv(os.path.join(DATADIR, f"base_mgate_preds_s{s}.csv"))
                    for s in SEEDS], ignore_index=True)
    bg = bg.groupby("name").agg(p_base=("p_base", "mean")).reset_index()
    return agg.merge(bg, on="name", how="left")


METHODS = [("SecProGNN", "p_base"), ("MKMP", "p_mkmp"), ("ESM-MLP", "p_esm"),
           ("加权融合", "p_weighted"), ("ReliFuse", "p_gated2")]
KS = [50, 100, 200, 500, 1000]


# ------------------------------------------------------------------ E18
def e18_planning(df):
    rows = []
    for K in KS:
        rec = {"候选数K": K, "肽段预算(2K)": 2 * K}
        for name, col in METHODS:
            top = df.sort_values(col, ascending=False).head(K)
            realized = top.y.mean()
            expected = top[col].mean()
            rec[f"{name} 预期精确率"] = round(float(expected), 4)
            rec[f"{name} 实际精确率"] = round(float(realized), 4)
            rec[f"{name} 规划误差"] = round(abs(float(expected) - float(realized)), 4)
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(EXP, "e18_planning_accuracy.csv"),
               index=False, encoding="utf-8-sig")
    return out


# ------------------------------------------------------------------ E19
def e19_optimal_depth(df, lambdas=(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)):
    """Net benefit (in units of confirmed detections) for verifying depth K.

    value per confirmed detection = 1; cost per peptide = lambda; a candidate
    costs 2 peptides -> net(K) = TP(K) - lambda * 2K.
    """
    grid = np.unique(np.concatenate([np.arange(5, 200, 5), np.arange(200, 1614, 20)]))
    rows = []
    for lam in lambdas:
        rec = {"单位验证成本 λ": lam}
        for name, col in METHODS:
            order = df.sort_values(col, ascending=False)
            tp = np.cumsum(order.y.values)[grid - 1]
            net = tp - lam * 2 * grid
            k_star = int(grid[int(np.argmax(net))])
            rec[f"{name} K*"] = k_star
            rec[f"{name} 净收益"] = round(float(net.max()), 1)
            rec[f"{name} K*处确认数"] = int(tp[int(np.argmax(net))])
        # oracle (uses labels): verify all true positives first
        order = df.sort_values("y", ascending=False)
        tp = np.cumsum(order.y.values)[grid - 1]
        net = tp - lam * 2 * grid
        rec["oracle K*"] = int(grid[int(np.argmax(net))])
        rec["oracle 净收益"] = round(float(net.max()), 1)
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(EXP, "e19_budget_optimality.csv"),
               index=False, encoding="utf-8-sig")
    return out, grid


# ------------------------------------------------------------------ fig 11
def fig11(df, plan, opt, grid):
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.9))
    cols = [C["grey"], C["blue"], C["green"], C["sky"], C["orange"]]
    names = ["SecProGNN", "MKMP", "ESM-MLP", "加权融合", "ReliFuse"]
    disp = ["SecProGNN", "MKMP", "ESM-MLP", "weighted", "ReliFuse"]

    # (a) planning error |expected - realised| at each budget
    ax = axes[0]
    x = np.arange(len(KS)); w = 0.16
    for j, (nm, dp, c) in enumerate(zip(names, disp, cols)):
        vals = [plan[f"{nm} 规划误差"].iloc[i] for i in range(len(KS))]
        ax.bar(x + (j - 2) * w, vals, w, color=c, edgecolor="white", lw=0.4, label=dp)
    ax.set_xticks(x); ax.set_xticklabels([f"Top-{k}" for k in KS], fontsize=6.4)
    ax.set(xlabel="Candidate budget", ylabel="|predicted − realised| precision",
           title="(a) Yield-planning error")
    ax.legend(fontsize=5.6, ncol=2)

    # (b) predicted vs realised precision
    ax2 = axes[1]
    for nm, dp, c, mk in zip(names, disp, cols, "osd^v"):
        xs = [plan[f"{nm} 预期精确率"].iloc[i] for i in range(len(KS))]
        ys = [plan[f"{nm} 实际精确率"].iloc[i] for i in range(len(KS))]
        ax2.scatter(xs, ys, s=16, color=c, marker=mk, label=dp, zorder=3)
    ax2.plot([0.45, 1.0], [0.45, 1.0], color="#BBBBBB", lw=0.9, ls="--")
    ax2.set(xlabel="Predicted precision of the shortlist",
            ylabel="Realised precision",
            title="(b) Calibrated yield forecasting")
    ax2.legend(fontsize=5.6)
    ax2.set_xlim(0.45, 1.02); ax2.set_ylim(0.45, 1.02)

    # (c) optimal verification depth vs cost ratio
    ax3 = axes[2]
    lam_axis = opt["单位验证成本 λ"].values
    for lab, col, c, ls in [("ReliFuse", "ReliFuse", C["orange"], "-"),
                            ("weighted", "加权融合", C["green"], "--"),
                            ("MKMP", "MKMP", C["blue"], ":"),
                            ("oracle", "oracle", C["grey"], "-.")]:
        ax3.plot(lam_axis, opt[f"{col} K*"].values, marker="o", ms=3.2, lw=1.4,
                 color=c, ls=ls, label=lab)
    ax3.set_yscale("log")
    ax3.set(xlabel="Verification cost per peptide  λ", ylabel="Recommended depth  K*",
            title="(c) Recommended verification depth")
    ax3.legend(fontsize=6.0)
    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig11_decision_policy")


if __name__ == "__main__":
    df = load()
    plan = e18_planning(df)
    print("=== E18 planning accuracy ===")
    print(plan[[c for c in plan.columns if c in ("候选数K",)] +
               [c for c in plan.columns if "规划误差" in c]].to_string(index=False))
    opt, grid = e19_optimal_depth(df)
    print("\n=== E19 optimal verification depth ===")
    print(opt[[c for c in opt.columns if c in ("单位验证成本 λ",) or "K*" in c]].to_string(index=False))
    fig11(df, plan, opt, grid)
