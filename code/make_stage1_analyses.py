# -*- coding: utf-8 -*-
"""Stage-1 analyses for v46.

E15  budget -> confirmed detections (peptide budget utility)
E16  strict homology-free evaluation (max identity to training set < 0.20 / <= 0.25)
E17  calibration head-to-head (weighted vs stacking baselines vs ReliFuse)
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps
import matplotlib.pyplot as plt
from sklearn.metrics import (roc_auc_score, matthews_corrcoef,
                             brier_score_loss, average_precision_score)

from style import C, save, DATADIR

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

SEEDS = [3407, 42, 123, 2024, 777, 1, 2, 3, 4, 5, 6, 7, 8, 10, 21, 99,
         202, 555, 888, 1000]
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, "..", "results"))
RNG = np.random.default_rng(20261002)


def ece(y, p, n_bins=10):
    y, p = np.asarray(y, float), np.asarray(p, float)
    edges = np.linspace(0, 1, n_bins + 1)
    e = 0.0
    for i in range(n_bins):
        m = ((p >= edges[i]) & (p < edges[i + 1])) if i < n_bins - 1 \
            else ((p >= edges[i]) & (p <= edges[i + 1]))
        if m.sum():
            e += m.mean() * abs(y[m].mean() - p[m].mean())
    return e


def load_main_preds():
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


# ------------------------------------------------------------------ E15
def e15_budget(df):
    methods = [("SecProGNN", "p_base"), ("MKMP", "p_mkmp"),
               ("加权融合", "p_weighted"), ("ReliFuse", "p_gated2")]
    rows = []
    for K in (50, 100, 200, 500, 1000):
        rec = {"候选数K": K, "肽段预算(2K)": 2 * K}
        tp = {}
        for name, col in methods:
            top = df.sort_values(col, ascending=False).head(K)
            tp[name] = int(top.y.sum())
            rec[f"{name} TP"] = tp[name]
            rec[f"{name} 精确率"] = round(tp[name] / K, 4)
            rec[f"{name} 每确认蛋白肽段数"] = round(2 * K / max(tp[name], 1), 2)
        rec["ΔTP ReliFuse−SecProGNN"] = tp["ReliFuse"] - tp["SecProGNN"]
        rec["ΔTP ReliFuse−MKMP"] = tp["ReliFuse"] - tp["MKMP"]
        rec["ΔTP ReliFuse−加权"] = tp["ReliFuse"] - tp["加权融合"]
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(EXP, "e15_budget_utility.csv"),
               index=False, encoding="utf-8-sig")
    return out


# ------------------------------------------------------------------ E16
def _metrics(y, p):
    return dict(AUC=roc_auc_score(y, p),
                PR_AUC=average_precision_score(y, p),
                MCC=matthews_corrcoef(y, (p >= .5).astype(int)),
                Brier=brier_score_loss(y, p), ECE=ece(y, p))


def e16_homology_free(df):
    # local fixed-seed RNG so this function is reproducible regardless of
    # any other RNG draws elsewhere in this script (v51 traceability fix)
    rng = np.random.default_rng(20261022)
    hi = pd.read_csv(os.path.join(DATADIR, "homology_identity.csv"))
    hi = hi.merge(df[["name", "p_esm", "p_weighted"]], on="name", how="left")
    subsets = [("全部测试集", hi.index == hi.index),
               ("同源一致度<0.20 (无同源)", hi.max_identity < 0.20),
               ("≤0.25", hi.max_identity <= 0.25),
               (">0.25 (含同源)", hi.max_identity > 0.25),
               (">0.40 (强同源)", hi.max_identity > 0.40)]
    methods = [("SecProGNN", "p_base"), ("MKMP", "p_mkmp"),
               ("ESM-MLP", "p_esm"), ("加权融合", "p_weighted"),
               ("ReliFuse", "p_lgf")]
    rows = []
    for lab, mask in subsets:
        sub = hi[np.asarray(mask)]
        y = sub.y.values
        rec = {"子集": lab, "n": len(sub), "阳性数": int(y.sum())}
        mcc = {}
        for name, col in methods:
            m = _metrics(y, sub[col].values)
            rec[f"{name} MCC"] = round(m["MCC"], 4)
            rec[f"{name} AUC"] = round(m["AUC"], 4)
            mcc[name] = m["MCC"]
        best = max(("SecProGNN", "MKMP", "ESM-MLP"), key=lambda k: mcc[k])
        rec["该子集最佳单路径"] = best
        rec["ΔMCC(ReliFuse−最佳单路径)"] = round(mcc["ReliFuse"] - mcc[best], 4)
        rec["ΔMCC(ReliFuse−SecProGNN)"] = round(mcc["ReliFuse"] - mcc["SecProGNN"], 4)
        # paired bootstrap CIs of ReliFuse - best single and ReliFuse - SecProGNN
        pb = sub[{"SecProGNN": "p_base", "MKMP": "p_mkmp",
                  "ESM-MLP": "p_esm"}[best]].values
        pc = sub.p_lgf.values
        ps = sub.p_base.values
        deltas, deltas_sec = [], []
        n = len(y)
        for _ in range(10000):
            idx = rng.integers(0, n, n)
            if len(np.unique(y[idx])) < 2:
                continue
            deltas.append(matthews_corrcoef(y[idx], (pc[idx] >= .5).astype(int))
                          - matthews_corrcoef(y[idx], (pb[idx] >= .5).astype(int)))
            deltas_sec.append(matthews_corrcoef(y[idx], (pc[idx] >= .5).astype(int))
                              - matthews_corrcoef(y[idx], (ps[idx] >= .5).astype(int)))
        rec["ΔMCC 95% CI"] = f"[{np.percentile(deltas,2.5):+.4f},{np.percentile(deltas,97.5):+.4f}]"
        rec["ΔMCC 95% CI(vs SecProGNN)"] = f"[{np.percentile(deltas_sec,2.5):+.4f},{np.percentile(deltas_sec,97.5):+.4f}]"
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(EXP, "e16_homology_free.csv"),
               index=False, encoding="utf-8-sig")
    return out


# ------------------------------------------------------------------ E17
def e17_calibration():
    res = pd.read_csv(os.path.join(DATADIR, "esm150m_v4_fusion_results.csv"))
    order = [("weighted", "加权融合"), ("stack", "logistic stacking"),
             ("stack_iso", "isotonic stacking"), ("gated2", "ReliFuse")]
    metrics = ["Brier", "ECE", "AUC", "MCC", "top500_tp"]
    rows = []
    for key, lab in order:
        sub = res[res.method == key]
        rows.append({"方法": lab, "Brier": sub.Brier.mean(), "ECE": sub.ECE.mean(),
                     "AUC": sub.AUC.mean(), "MCC": sub.MCC.mean(),
                     "Top500": sub.top500_tp.mean()})
    tbl = pd.DataFrame(rows).round(4)
    tbl.to_csv(os.path.join(EXP, "e17_calibration_headtohead.csv"),
               index=False, encoding="utf-8-sig")
    # paired tests: ReliFuse vs each baseline
    test_rows = []
    base = res[res.method == "gated2"]
    for key, lab in order[:-1]:
        other = res[res.method == key]
        pvals = []
        for met in metrics:
            a = base.sort_values("seed")[met].values
            b = other.sort_values("seed")[met].values
            try:
                _, p = sps.wilcoxon(a, b)
            except Exception:
                p = np.nan
            pvals.append(p)
        order_idx = np.argsort(pvals)
        adj = np.empty(len(metrics))
        m = len(metrics)
        for rank, i in enumerate(order_idx):
            adj[i] = min(1.0, (m - rank) * pvals[i])
        for rank, i in enumerate(order_idx[1:], 1):
            adj[i] = max(adj[i], adj[order_idx[rank - 1]])
        for met, p, pa in zip(metrics, pvals, adj):
            a = base.sort_values("seed")[met].values
            b = other.sort_values("seed")[met].values
            d = a - b
            better = "lower" if met in ("Brier", "ECE") else "higher"
            win = int((d < 0).sum()) if better == "lower" else int((d > 0).sum())
            lose = int((d > 0).sum()) if better == "lower" else int((d < 0).sum())
            test_rows.append({"对比": f"ReliFuse vs {lab}", "指标": met,
                              "均值A": round(a.mean(), 4), "均值B": round(b.mean(), 4),
                              "Δ": round(d.mean(), 4), "胜/负": f"{win}/{lose}",
                              "Wilcoxon p": f"{p:.3g}", "Holm p": f"{pa:.3g}"})
    tt = pd.DataFrame(test_rows)
    tt.to_csv(os.path.join(EXP, "e17_calibration_paired.csv"),
              index=False, encoding="utf-8-sig")
    return tbl, tt


# ------------------------------------------------------------------ figures
def fig8_v2(df):
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.8))
    y = df.y.values
    ax = axes[0]
    for lab, col, c in [("MKMP", "p_mkmp", C["blue"]),
                        ("ReliFuse", "p_gated2", C["orange"])]:
        p = df[col].values
        edges = np.linspace(0, 1, 11)
        xs, ys = [], []
        for i in range(10):
            m = (p >= edges[i]) & (p < edges[i + 1]) if i < 9 else (p >= edges[i]) & (p <= edges[i + 1])
            if m.sum():
                xs.append(p[m].mean()); ys.append(y[m].mean())
        ax.plot(xs, ys, marker="o", ms=3, lw=1.5, color=c,
                label=f"{lab} (ECE {ece(y, p):.3f})")
    ax.plot([0, 1], [0, 1], color="#BBBBBB", lw=0.8, ls="--")
    ax.set(xlabel="Predicted probability", ylabel="Observed frequency",
           title="(a) Reliability")
    ax.legend(fontsize=6.0, loc="upper left")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)

    ax2 = axes[1]
    dca = pd.read_csv(os.path.join(DATADIR, "dca_netbenefit.csv"))
    for lab, c, mk in [("MKMP", C["blue"], "o"), ("weighted", C["green"], "s"),
                       ("ReliFuse", C["orange"], "^")]:
        ax2.plot(dca.threshold, dca[lab], marker=mk, ms=3, lw=1.5, color=c, label=lab)
    ax2.axhline(0, color="#BBBBBB", lw=0.8)
    ax2.set(xlabel="Threshold probability", ylabel="Net benefit",
            title="(b) Decision curve analysis")
    ax2.legend(fontsize=6.2)

    budget = e15_budget(df)
    ax3 = axes[2]
    K = budget["候选数K"].values
    pep = budget["肽段预算(2K)"].values
    for name, c, mk in [("SecProGNN", C["grey"], "d"), ("MKMP", C["blue"], "o"),
                        ("ReliFuse", C["orange"], "^")]:
        ax3.plot(pep, budget[f"{name} TP"].values, marker=mk, ms=3.5, lw=1.6,
                 color=c, ls="--" if name == "SecProGNN" else "-", label=name)
    for xi, yi, d in zip(pep, budget["ReliFuse TP"].values,
                         budget["ΔTP ReliFuse−SecProGNN"].values):
        ax3.annotate(f"+{int(d)}", (xi, yi), textcoords="offset points",
                     xytext=(0, 7), ha="center", fontsize=6.0, color=C["orange"])
    ax3.set_xscale("log")
    ax3.set_xticks(pep); ax3.set_xticklabels([str(p) for p in pep])
    ax3.set(xlabel="Peptide budget (2 peptides / candidate)",
            ylabel="Confirmed detections", title="(c) Budget → confirmed detections")
    ax3.legend(fontsize=6.0, loc="upper left")

    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig8_calibration_decision")


def fig9_v2(df, hom):
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9))
    cohorts = ["Main", "GlycoHBF", "Cohort A", "Cohort B"]
    mkmp = [0.9588, 0.8215, 0.513, 0.5896]
    esm = [0.9548, 0, 0.845, 0.6684]
    clf = [0.9687, 0.8718, 0.520, 0.5908]
    x = np.arange(4); w = 0.27
    ax = axes[0]
    ax.bar(x - w, mkmp, w, color=C["blue"], edgecolor="white", lw=0.5, label="MKMP")
    ax.bar(x, esm, w, color=C["green"], edgecolor="white", lw=0.5, label="ESM-MLP")
    ax.bar(x + w, clf, w, color=C["orange"], edgecolor="white", lw=0.5, label="ReliFuse")
    ax.set_xticks(x); ax.set_xticklabels(cohorts, fontsize=6.6, rotation=12)
    ax.set(ylabel="AUC", title="(a) Across cohorts and body fluids")
    ax.legend(fontsize=5.8, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.24))
    ax.set_ylim(0, 1.05)

    ax2 = axes[1]
    hi = pd.read_csv(os.path.join(DATADIR, "homology_identity.csv"))
    edges = [0, 0.2, 0.4, 0.6, 0.8, 1.01]
    labs = ["<0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", ">0.8"]
    dm, lo, hh, ns = [], [], [], []
    for i in range(5):
        m = (hi.max_identity >= edges[i]) & (hi.max_identity < edges[i + 1])
        sub = hi[m]
        ns.append(len(sub))
        if sub.y.nunique() < 2:
            dm.append(0.0); lo.append(0.0); hh.append(0.0); continue
        a = matthews_corrcoef(sub.y, (sub.p_lgf >= .5).astype(int))
        b = matthews_corrcoef(sub.y, (sub.p_mkmp >= .5).astype(int))
        dm.append(a - b)
        n = len(sub); boot = []
        for _ in range(1000):
            idx = RNG.integers(0, n, n)
            if len(np.unique(sub.y.values[idx])) < 2:
                continue
            boot.append(matthews_corrcoef(sub.y.values[idx], (sub.p_lgf.values[idx] >= .5).astype(int))
                        - matthews_corrcoef(sub.y.values[idx], (sub.p_mkmp.values[idx] >= .5).astype(int)))
        lo.append(np.percentile(boot, 2.5)); hh.append(np.percentile(boot, 97.5))
    err = [np.array(dm) - np.array(lo), np.array(hh) - np.array(dm)]
    bars = ax2.bar(labs, dm, yerr=err, capsize=3,
                   color=[C["grey"]] * 4 + [C["orange"]], edgecolor=C["grey"], lw=0.5)
    ax2.axhline(0, color="#BBBBBB", lw=0.8)
    ax2.set(xlabel="Max. sequence identity to training set",
            ylabel=r"$\Delta$MCC (fusion − structure)",
            title="(b) Fusion benefit vs homology")
    ax2.bar_label(bars, labels=[f"n={n}" for n in ns], fontsize=5.8, padding=1)
    ax2.tick_params(axis="x", labelsize=6.4)

    ax3 = axes[2]
    cc = pd.read_csv(os.path.join(DATADIR, "clustercv_esm150m_fusion_results.csv"))
    order = ["mkmp", "esm", "weighted", "gated2"]
    labs2 = ["MKMP", "ESM", "weighted", "ReliFuse"]
    cols = [C["blue"], C["green"], C["grey"], C["orange"]]
    mu = [cc[cc.method == o].MCC.mean() for o in order]
    sd = [cc[cc.method == o].MCC.std(ddof=1) for o in order]
    ax3.bar(labs2, mu, yerr=sd, capsize=3, color=cols, edgecolor="white", lw=0.5)
    ax3.set(ylabel="MCC (cluster-disjoint CV)",
            title="(c) Homology-aware cross-validation")
    ax3.tick_params(axis="x", labelsize=6.0)
    ax3.set_xticklabels(labs2, rotation=15)
    ax3.set_ylim(0.70, 0.85)
    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig9_transfer_robustness")


if __name__ == "__main__":
    df = load_main_preds()
    b = e15_budget(df)
    print("=== E15 budget utility ===")
    print(b.to_string(index=False))
    h = e16_homology_free(df)
    print("\n=== E16 homology-free subsets ===")
    print(h.to_string(index=False))
    tbl, tt = e17_calibration()
    print("\n=== E17 calibration head-to-head ===")
    print(tbl.to_string(index=False))
    print(tt.to_string(index=False))
    fig8_v2(df)
    hi = pd.read_csv(os.path.join(DATADIR, "homology_identity.csv"))
    fig9_v2(df, hi)
