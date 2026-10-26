# -*- coding: utf-8 -*-
"""Data-driven v46 figures (F3-F9 and supplementary panels)."""

import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             matthews_corrcoef, f1_score, brier_score_loss,
                             precision_recall_curve, roc_curve, calinski_harabasz_score)

from style import C, save, DATADIR, panel_label

SEEDS = [3407, 42, 123, 2024, 777, 1, 2, 3, 4, 5, 6, 7, 8, 10, 21, 99,
         202, 555, 888, 1000]
ESM_SEEDS = [3407, 42, 123, 2024, 777]
RNG = np.random.default_rng(20261002)


def ece(y, p, n_bins=10):
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    edges = np.linspace(0, 1, n_bins + 1)
    e = 0.0
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        m = (p >= lo) & (p < hi) if i < n_bins - 1 else (p >= lo) & (p <= hi)
        if m.sum():
            e += m.mean() * abs(y[m].mean() - p[m].mean())
    return e


def load_main():
    fr = [pd.read_csv(os.path.join(DATADIR, f"preds_esm150m_v4_s{s}.csv"))
          for s in SEEDS]
    ap = pd.concat(fr, ignore_index=True)
    agg = ap.groupby("name").agg(
        y=("y", "max"), p_mkmp=("p_mkmp", "mean"), p_esm=("p_esm", "mean"),
        p_gated2=("p_gated2", "mean"), p_weighted=("p_weighted", "mean"),
        p_gated=("p_gated", "mean"), z_mkmp=("z_mkmp_cal", "mean"),
        z_esm=("z_esm_cal", "mean"), g=("g_gated2", "mean"),
        b=("b_gated2", "mean")).reset_index()
    bg = pd.concat([pd.read_csv(os.path.join(DATADIR, f"base_mgate_preds_s{s}.csv"))
                    for s in SEEDS], ignore_index=True)
    bg = bg.groupby("name").agg(p_base=("p_base", "mean"),
                                p_mgate=("p_mgate", "mean")).reset_index()
    df = agg.merge(bg, on="name", how="left")
    pf = np.load(os.path.join(DATADIR, "plddt_features_8069.npy"))
    pn = np.load(os.path.join(DATADIR, "plddt_names_8069.npy"), allow_pickle=True)
    pmap = {str(n): float(pf[i, 0]) for i, n in enumerate(pn)}
    df["plddt"] = df["name"].map(pmap)
    return df.dropna(subset=["plddt"]).reset_index(drop=True)


# ------------------------------------------------------------------ fig 3
def fig3(df):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    series = [("SecProGNN", "p_base", C["grey"], (0, (3, 2))),
              ("MKMP (structure)", "p_mkmp", C["blue"], "-"),
              ("ESM-MLP (sequence)", "p_esm", C["green"], (0, (5, 1))),
              ("weighted fusion", "p_weighted", C["sky"], (0, (1, 1))),
              ("ReliFuse (fusion)", "p_gated2", C["orange"], "-")]
    y = df.y.values
    for lab, col, c, ls in series:
        p = df[col].values
        fpr, tpr, _ = roc_curve(y, p)
        pre, rec, _ = precision_recall_curve(y, p)
        axes[0].plot(fpr, tpr, color=c, ls=ls, lw=1.2,
                     label=f"{lab} (AUC {roc_auc_score(y,p):.3f})")
        axes[1].plot(rec, pre, color=c, ls=ls, lw=1.2,
                     label=f"{lab} (AP {average_precision_score(y,p):.3f})")
    axes[0].plot([0, 1], [0, 1], color="#BBBBBB", lw=0.8, ls=":")
    axes[0].set(xlabel="False positive rate", ylabel="True positive rate",
                title="(a) ROC — main benchmark")
    axes[1].set(xlabel="Recall", ylabel="Precision",
                title="(b) Precision–recall")
    axes[0].legend(loc="lower right", fontsize=5.6)
    axes[1].legend(loc="lower left", fontsize=5.6)
    axes[0].set_xlim(-0.01, 1.01); axes[0].set_ylim(-0.01, 1.01)
    fig.tight_layout()
    save(fig, "fig3_main_roc_pr")


# ------------------------------------------------------------------ fig 4
def fig4(df):
    d = df.sort_values("plddt").reset_index(drop=True)
    y = d.y.values
    dec = pd.qcut(d.plddt, 6, labels=False, duplicates="drop")
    xs2, curves = [], {k: [] for k in ("p_mkmp", "p_esm", "p_gated2")}
    dlt, lo, hi = [], [], []
    for b in sorted(pd.unique(dec)):
        m = (dec == b).values
        yy = y[m]
        xs2.append(d.plddt.values[m].mean())
        for k in curves:
            curves[k].append(matthews_corrcoef(yy, (d[k].values[m] >= .5).astype(int)))
        pg, pm, pe = (d.p_gated2.values[m], d.p_mkmp.values[m], d.p_esm.values[m])
        dlt.append(matthews_corrcoef(yy, (pg >= .5).astype(int))
                   - max(matthews_corrcoef(yy, (pm >= .5).astype(int)),
                         matthews_corrcoef(yy, (pe >= .5).astype(int))))
        boot = []
        n = len(yy)
        for _ in range(400):
            idx = RNG.integers(0, n, n)
            if len(np.unique(yy[idx])) < 2:
                continue
            boot.append(matthews_corrcoef(yy[idx], (pg[idx] >= .5).astype(int))
                        - max(matthews_corrcoef(yy[idx], (pm[idx] >= .5).astype(int)),
                              matthews_corrcoef(yy[idx], (pe[idx] >= .5).astype(int))))
        lo.append(np.percentile(boot, 2.5)); hi.append(np.percentile(boot, 97.5))

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))
    ax = axes[0]
    for (lab, k, c) in [("MKMP", "p_mkmp", C["blue"]), ("ESM-MLP", "p_esm", C["green"]),
                        ("ReliFuse", "p_gated2", C["orange"])]:
        ax.plot(xs2, curves[k], marker="o", ms=3, color=c, lw=1.6, label=lab)
    ax.set(xlabel="pLDDT (bin mean)", ylabel="MCC",
           title="(a) Performance vs pLDDT")
    ax.legend(loc="lower right")
    ax2 = axes[1]
    ax2.axhline(0, color="#BBBBBB", lw=0.8)
    ax2.fill_between(xs2, lo, hi, color=C["orange"], alpha=0.20, lw=0)
    ax2.plot(xs2, dlt, marker="o", ms=3, color=C["orange"], lw=1.8)
    ax2.axhline(0, color="#BBBBBB", lw=0.8)
    ax2.set(xlabel="pLDDT (bin mean)", ylabel=r"$\Delta$MCC (fusion − best single)",
            title="(b) Fusion gain vs pLDDT")
    ax3 = axes[2]
    dis = np.abs(df.p_mkmp - df.p_esm)
    q = pd.qcut(dis, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    yy = df.y.values
    pf = (df.p_gated2.values >= .5).astype(int)
    pm = (df.p_mkmp.values >= .5).astype(int)
    pe = (df.p_esm.values >= .5).astype(int)
    gains = []
    for qq in ["Q1", "Q2", "Q3", "Q4"]:
        m = (q == qq).values
        gains.append((pf[m] == yy[m]).mean() - max((pm[m] == yy[m]).mean(), (pe[m] == yy[m]).mean()))
    bars = ax3.bar(["Q1", "Q2", "Q3", "Q4"], gains, color=[C["lightgrey"]] * 3 + [C["orange"]],
                   edgecolor=C["grey"], lw=0.6)
    ax3.set(xlabel="Path-disagreement quartile", ylabel="Accuracy gain over best single",
            title="(c) Gain by disagreement quartile")
    ax3.bar_label(bars, fmt="%+.3f", fontsize=6.2, padding=1)
    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig4_complementarity")


# ------------------------------------------------------------------ fig 5
def fig5(df):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
    ax = axes[0]
    hb = ax.hexbin(df.p_esm, df.p_mkmp, gridsize=42, cmap="Blues",
                   norm=LogNorm(vmin=1), mincnt=1, linewidths=0)
    ax.axvline(0.5, color=C["grey"], lw=0.8, ls="--")
    ax.axhline(0.5, color=C["grey"], lw=0.8, ls="--")
    ax.set(xlabel="Sequence-path probability $p_{ESM}$",
           ylabel="Structure-path probability $p_{MKMP}$",
           title="(a) Two-path probability map (20-seed aggregate)")
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("proteins", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    ax.text(0.06, 0.94, "agree", fontsize=6.5, color=C["grey"], transform=ax.transAxes)
    ax.text(0.62, 0.06, "agree", fontsize=6.5, color=C["grey"], transform=ax.transAxes)
    ax.text(0.30, 0.05, "structure only", fontsize=6.5, color=C["blue"],
            transform=ax.transAxes)
    ax.text(0.60, 0.88, "sequence only", fontsize=6.5, color=C["green"],
            transform=ax.transAxes)
    ax.grid(False)

    ax2 = axes[1]
    y = df.y.values
    pf = (df.p_gated2.values >= .5).astype(int)
    pm = (df.p_mkmp.values >= .5).astype(int)
    pe = (df.p_esm.values >= .5).astype(int)
    labels, ns, gains = [], [], []
    for a, b, lab in [(0, 0, "both 0"), (1, 1, "both 1"),
                      (0, 1, "ESM only"), (1, 0, "MKMP only")]:
        m = (pm == a) & (pe == b)
        labels.append(lab); ns.append(int(m.sum()))
        gains.append((pf[m] == y[m]).mean() - max((pm[m] == y[m]).mean(), (pe[m] == y[m]).mean()))
    cols = [C["lightgrey"], C["lightgrey"], C["green"], C["blue"]]
    bars = ax2.bar(labels, gains, color=cols, edgecolor=C["grey"], lw=0.6)
    ax2.axhline(0, color="#BBBBBB", lw=0.8)
    ax2.set(ylabel="Fusion accuracy − best single path",
            title="(b) Fusion benefit by decision quadrant")
    ax2.bar_label(bars, labels=[f"{g:+.3f}\nn={n}" for g, n in zip(gains, ns)],
                  fontsize=6.2, padding=2)
    ax2.set_ylim(min(-0.02, min(gains) - 0.01), max(gains) * 1.40)
    ax2.tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    save(fig, "fig5_decision_map")


# ------------------------------------------------------------------ fig 6
def fig6():
    trs = []
    for s in ESM_SEEDS:
        t = pd.read_csv(os.path.join(DATADIR, f"trajectory_v12_s{s}.csv"))
        keep = ["name"] + [f"a{l}_{k}" for l in range(3) for k in range(4)]
        trs.append(t[keep])
    tr = pd.concat(trs).groupby("name").mean().reset_index()
    df = load_main()[["name", "y", "plddt"]].merge(tr, on="name", how="inner")
    df["stratum"] = pd.cut(df.plddt, [0, 70, 90, 100],
                           labels=["pLDDT < 70", "70–90", "> 90"])
    kmap = {0: r"$\kappa_1$ cosine", 1: r"$\kappa_2$ chain",
            2: r"$\kappa_3$ degree", 3: r"$\kappa_4$ RBF"}
    kcols = [C["sky"], C["green"], C["amber"], C["blue"]]

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.8))
    ax = axes[0]
    strata = ["pLDDT < 70", "70–90", "> 90"]
    x = np.arange(len(strata)); w = 0.2
    for k in range(4):
        vals, errs = [], []
        for st in strata:
            v = df[df.stratum == st][[f"a{l}_{k}" for l in range(3)]].values.mean(1)
            vals.append(v.mean())
            errs.append(1.96 * v.std(ddof=1) / np.sqrt(len(v)))
        ax.bar(x + (k - 1.5) * w, vals, w, yerr=errs, capsize=1.6,
               color=kcols[k], edgecolor="white", lw=0.4, label=kmap[k])
    ax.set_xticks(x); ax.set_xticklabels(["< 70", "70–90", "> 90"])
    ax.set(xlabel="pLDDT stratum", ylabel="Mean kernel weight  $\\alpha_k$",
           title="(a) Kernel usage by structural confidence")
    ax.legend(fontsize=5.8, ncol=2, loc="upper center")
    ax.set_ylim(0.15, 0.33)

    ax2 = axes[1]
    for k in range(4):
        ys = [df[f"a{l}_{k}"].mean() for l in range(3)]
        es = [1.96 * df[f"a{l}_{k}"].std(ddof=1) / np.sqrt(len(df)) for l in range(3)]
        ax2.errorbar([0, 1, 2], ys, yerr=es, marker="o", ms=3.5, lw=1.5,
                     capsize=2, color=kcols[k], label=kmap[k])
    ax2.set_xticks([0, 1, 2]); ax2.set_xticklabels(["layer 1", "layer 2", "layer 3"])
    ax2.set(xlabel="Message-passing depth", ylabel="Mean kernel weight  $\\alpha_k$",
            title="(b) Kernel weighting across layers")
    ax2.legend(fontsize=5.8, ncol=2)

    ax3 = axes[2]
    bc = pd.read_csv(os.path.join(DATADIR, "mkmp_v12_bio_corr.csv"))
    bc = bc[bc.feature == "plddt"].sort_values("mech")
    r = bc.r.values
    bars = ax3.barh([kmap[int(m)] for m in bc.mech.values], r,
                    color=[kcols[int(m)] for m in bc.mech.values],
                    edgecolor=C["grey"], lw=0.5)
    ax3.axvline(0, color="#BBBBBB", lw=0.8)
    ax3.set(xlabel="Pearson $r$ ($\\alpha_k$ vs pLDDT)",
            title="(c) Kernel weight vs structural confidence")
    ax3.bar_label(bars, fmt="%+.2f", fontsize=6.2, padding=2)
    ax3.set_xlim(-0.62, 0.62)
    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig6_mechanism_kernels")


# ------------------------------------------------------------------ fig 7
def fig7(df):
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.8))
    ax = axes[0]
    hb = ax.hexbin(df.g, np.clip(df.b, -4, 4), gridsize=40, cmap="Oranges",
                   norm=LogNorm(vmin=1), mincnt=1, linewidths=0)
    ax.set(xlabel="Fusion weight  $g$", ylabel="Per-protein bias  $b$  (clipped to ±4)",
           title="(a) Fusion-operator parameters")
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=6); cb.set_label("proteins", fontsize=6.5)
    ax.grid(False)

    ax2 = axes[1]
    contrib_g = np.abs(df.g.values * (df.z_esm.values - df.z_mkmp.values))
    contrib_b = np.abs(df.b.values)
    bins = np.linspace(0, 10, 40)
    ax2.hist(contrib_g, bins=bins, histtype="step", lw=1.6, color=C["orange"],
             label=f"$|g\\,(z_{{2c}}-z_{{1c}})|$  mean {contrib_g.mean():.2f}")
    ax2.hist(contrib_b, bins=bins, histtype="step", lw=1.6, color=C["blue"],
             label=f"$|b|$  mean {contrib_b.mean():.2f}")
    ax2.set(xlabel="Magnitude of term in $z_{fuse}$", ylabel="Proteins",
            title="(b) Weight vs bias contribution")
    ax2.legend(fontsize=6.0)

    ax3 = axes[2]
    cf = pd.read_csv(os.path.join(DATADIR, "e10_counterfactual.csv"))
    pert = ["pLDDT=30", "pLDDT=50", "pLDDT=70", "pLDDT=90"]
    g = [cf[cf.perturbation == p].mean_dg.values for p in pert]
    b = [cf[cf.perturbation == p].mean_db.values for p in pert]
    x = np.arange(4); w = 0.36
    ax3.bar(x - w / 2, [v.mean() for v in g], w,
            yerr=[1.96 * v.std(ddof=1) / np.sqrt(len(v)) for v in g],
            capsize=2, color=C["orange"], label=r"$\Delta g$")
    ax3.bar(x + w / 2, [v.mean() for v in b], w,
            yerr=[1.96 * v.std(ddof=1) / np.sqrt(len(v)) for v in b],
            capsize=2, color=C["blue"], label=r"$\Delta b$")
    ax3.axhline(0, color="#BBBBBB", lw=0.8)
    ax3.set_xticks(x); ax3.set_xticklabels(["30", "50", "70", "90"])
    ax3.set(xlabel="Injected pLDDT value", ylabel="Mean response ($\\Delta$ logit)",
            title="(c) Counterfactual response")
    ax3.legend(fontsize=6.5)
    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig7_fusion_anatomy")


# ------------------------------------------------------------------ fig 8
def fig8(df):
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))
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
    ax2.legend(fontsize=6.5)

    ax3 = axes[2]
    tk = pd.read_csv(os.path.join("/".join(DATADIR.split("\\")[:-1]), "results",
                                  "e5_topk_budget.csv"))
    ks = [50, 100, 200, 500, 1000]
    ax3.plot(ks, [tk["Top50 TP"][0], tk["Top100 TP"][0], tk["Top200 TP"][0],
                  tk["Top500 TP"][0], tk["Top1000 TP"][0]], marker="o", ms=3.5,
             lw=1.6, color=C["blue"], label="MKMP")
    ax3.plot(ks, [tk["Top50 TP"][1], tk["Top100 TP"][1], tk["Top200 TP"][1],
                  tk["Top500 TP"][1], tk["Top1000 TP"][1]], marker="^", ms=3.5,
             lw=1.6, color=C["orange"], label="ReliFuse")
    ax3.plot(ks, [tk["Top50 TP"][2], tk["Top100 TP"][2], tk["Top200 TP"][2],
                  tk["Top500 TP"][2], tk["Top1000 TP"][2]], marker="s", ms=3.5,
             lw=1.2, color=C["green"], ls="--", label="weighted")
    ax3.set_xscale("log")
    ax3.set_xticks(ks); ax3.set_xticklabels([str(k) for k in ks])
    ax3.set(xlabel="Candidate budget $K$", ylabel="Cumulative true positives",
            title="(c) Ranking under fixed budget")
    ax3.legend(fontsize=6.5, loc="lower right")
    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig8_calibration_decision")


# ------------------------------------------------------------------ fig 9
def fig9(df):
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9))
    # (a) transfer AUC
    cohorts = ["Main", "GlycoHBF", "Cohort A", "Cohort B"]
    mkmp = [0.9588, 0.8215, 0.513, 0.5896]
    esm = [0.9548, None, 0.845, 0.6684]
    clf = [0.9687, 0.8718, 0.520, 0.5908]
    x = np.arange(4); w = 0.27
    ax = axes[0]
    ax.bar(x - w, mkmp, w, color=C["blue"], edgecolor="white", lw=0.5, label="MKMP")
    ax.bar(x, [e if e is not None else 0 for e in esm], w, color=C["green"],
           edgecolor="white", lw=0.5, label="ESM-MLP")
    ax.bar(x + w, clf, w, color=C["orange"], edgecolor="white", lw=0.5, label="ReliFuse")
    ax.set_xticks(x); ax.set_xticklabels(cohorts, fontsize=6.6, rotation=12)
    ax.set(ylabel="AUC", title="(a) Across cohorts and body fluids")
    ax.legend(fontsize=6.0, loc="upper center", ncol=3,
              bbox_to_anchor=(0.5, -0.22))
    ax.set_ylim(0, 1.05)

    # (b) homology
    ax2 = axes[1]
    hi = pd.read_csv(os.path.join(DATADIR, "homology_identity.csv"))
    hi = hi.dropna(subset=["max_identity"])
    edges = [0, 0.2, 0.4, 0.6, 0.8, 1.01]
    labs = ["<0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", ">0.8"]
    dm, ns = [], []
    for i in range(5):
        m = (hi.max_identity >= edges[i]) & (hi.max_identity < edges[i + 1])
        sub = hi[m]
        if sub.y.nunique() < 2:
            dm.append(np.nan); ns.append(len(sub)); continue
        a = matthews_corrcoef(sub.y, (sub.p_lgf >= .5).astype(int))
        b = matthews_corrcoef(sub.y, (sub.p_mkmp >= .5).astype(int))
        dm.append(a - b); ns.append(len(sub))
    bars = ax2.bar(labs, dm, color=[C["grey"]] * 4 + [C["orange"]],
                   edgecolor=C["grey"], lw=0.5)
    ax2.axhline(0, color="#BBBBBB", lw=0.8)
    ax2.set(xlabel="Max. sequence identity to training set",
            ylabel=r"$\Delta$MCC (fusion − structure)",
            title="(b) Fusion benefit vs homology")
    ax2.bar_label(bars, labels=[f"n={n}" for n in ns], fontsize=5.8, padding=1)
    ax2.tick_params(axis="x", labelsize=6.4)

    # (c) cluster CV
    ax3 = axes[2]
    cc = pd.read_csv(os.path.join(DATADIR, "clustercv_esm150m_fusion_results.csv"))
    order = ["mkmp", "esm", "weighted", "gated2"]
    labs2 = ["MKMP", "ESM-MLP", "weighted", "ReliFuse"]
    cols = [C["blue"], C["green"], C["grey"], C["orange"]]
    mu = [cc[cc.method == o].MCC.mean() for o in order]
    sd = [cc[cc.method == o].MCC.std(ddof=1) for o in order]
    ax3.bar(labs2, mu, yerr=sd, capsize=3, color=cols, edgecolor="white", lw=0.5)
    ax3.set(ylabel="MCC (cluster-disjoint CV)",
            title="(c) Homology-aware cross-validation")
    ax3.tick_params(axis="x", labelsize=6.0)
    ax3.set_xticklabels(["MKMP", "ESM", "weighted", "ReliFuse"], rotation=15)
    ax3.set_ylim(0.70, 0.85)
    for a in axes:
        a.grid(alpha=0.18)
    fig.tight_layout()
    save(fig, "fig9_transfer_robustness")


# ---------------------------------------------------------- supplementary
def supp_gate():
    e12 = pd.read_csv(os.path.join(DATADIR, "e12_gate_identifiability.csv"))
    nice = {"learned": "learned\n$g,b$", "fixed_gbar": "fixed\n$g=\\bar g$",
            "fixed_g05_b": "fixed $g$=0.5\n(keep $b$)", "fixed_g05_nob": "fixed $g$=0.5\n(no $b$)",
            "perm_g": "permuted\n$g$", "perm_gb": "permuted\n$g$ and $b$"}
    e12["lab"] = e12.variant.map(lambda v: nice.get(v, v))
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    cols = [C["orange"]] + [C["grey"]] * (len(e12) - 1)
    axes[0].bar(e12.lab, e12.MCC_mean, yerr=e12.MCC_sd, capsize=2.5,
                color=cols, edgecolor="white", lw=0.5)
    axes[0].set(ylabel="MCC", title="(a) Gate identifiability — discrimination")
    axes[0].tick_params(axis="x", labelsize=5.8)
    axes[0].set_ylim(0.65, 0.83)
    axes[1].bar(e12.lab, e12.ECE_mean, yerr=e12.ECE_sd, capsize=2.5,
                color=cols, edgecolor="white", lw=0.5)
    axes[1].set(ylabel="ECE", title="(b) Gate identifiability — calibration")
    axes[1].tick_params(axis="x", labelsize=5.8)
    fig.tight_layout()
    save(fig, "supp_gate_identifiability")


def supp_ladder():
    res = pd.read_csv(os.path.join(DATADIR, "esm150m_v4_fusion_results.csv"))
    order = ["esm", "mkmp", "stack", "weighted", "gated", "gated2"]
    labs = ["ESM-MLP", "MKMP", "logistic\nstacking", "weighted", "gated", "ReliFuse"]
    m = res.groupby("method")[["AUC", "ECE", "MCC"]].mean().loc[order]
    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    xs = np.arange(len(order))
    ax.plot(xs, m.AUC, marker="o", ms=4, lw=1.7, color=C["blue"], label="AUC")
    ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=6.2)
    ax.set_ylabel("AUC", color=C["blue"])
    ax.set_ylim(0.945, 0.975)
    ax2 = ax.twinx()
    ax2.plot(xs, m.ECE, marker="s", ms=3.6, lw=1.5, ls="--", color=C["orange"], label="ECE")
    ax2.set_ylabel("ECE", color=C["orange"])
    ax2.set_ylim(0.015, 0.070)
    ax2.grid(False)
    ax2.spines["right"].set_visible(True)
    for i, (a, e) in enumerate(zip(m.AUC, m.ECE)):
        ax.annotate(f"{a:.4f}", (i, a), textcoords="offset points", xytext=(0, 6),
                    ha="center", fontsize=5.8, color=C["blue"])
        ax2.annotate(f"{e:.3f}", (i, e), textcoords="offset points", xytext=(0, -10),
                     ha="center", fontsize=5.8, color=C["orange"])
    ax.set_title("Pipeline stage ladder")
    fig.tight_layout()
    save(fig, "supp_pipeline_ladder")


if __name__ == "__main__":
    df = load_main()
    print("main benchmark proteins:", len(df))
    fig3(df)
    fig4(df)
    fig5(df)
    fig6()
    fig7(df)
    fig8(df)
    fig9(df)
    supp_gate()
    supp_ladder()
