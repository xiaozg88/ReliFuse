# -*- coding: utf-8 -*-
"""Regenerate the v45 figures that carried stale naming / stale numbers.

Produces into paper_v39_v40/v45_figures/:
  fig4_plddt_threeway.png      corrected pLDDT strata (n=531/848/234)
  fig_go_enrichment.png        "MKMP-LGF" -> "ReliFuse"
  fig_bias_distribution.png    "MKMP-LGF" -> "ReliFuse"
  fig_e2_overlap_strata.png    legend "MKMP-LGF" -> "ReliFuse"
  fig_e8_label_sensitivity.png legend "MKMP-LGF" -> "ReliFuse"
"""

import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# Development-tree root, four levels above this file. This provenance script
# reads upstream v1 pipeline outputs that are not part of the release
# repository; released counterparts live under ../data/ and ../results/.
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
V1 = os.path.join(ROOT, "v1")
OUT = os.path.join(ROOT, "v2", "paper_v39_v40", "v48", "figures")
RES0924 = os.path.join(V1, "results", "0924_experiments")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial"],
    "axes.unicode_minus": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
})


def load_agg(prob_cols=("p_mkmp", "p_esm", "p_gated2")):
    files = sorted(glob.glob(os.path.join(
        V1, "results", "logit_fusion", "preds_esm150m_v4_s*.csv")))
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    cols = {"y": "max", "g_gated2": "mean", "b_gated2": "mean"}
    for c in prob_cols:
        cols[c] = "mean"
    return df.groupby("name").agg(cols).reset_index()


# ---------------------------------------------------------------- fig 4
def fig4():
    strata = ["pLDDT < 70\n(n=531, 6.0% pos)",
              "70 \u2264 pLDDT \u2264 90\n(n=848, 33.7% pos)",
              "pLDDT > 90\n(n=234, 70.1% pos)"]
    mkmp = [0.6026, 0.7234, 0.7816]
    esm = [0.6536, 0.7371, 0.7146]
    clf = [0.7901, 0.7538, 0.8135]
    # Delta vs the BEST single path in each stratum (v45 corrected definition).
    deltas = ["+0.137\nCI [.019,.214]",
              "+0.017\nCI [\u2212.022,.042]",
              "+0.032\nCI [\u2212.028,.086]"]

    x = np.arange(3)
    w = 0.24
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    b1 = ax.bar(x - w, mkmp, w, label="MKMP (structure)", color="#1f77b4",
                edgecolor="white", zorder=3)
    b2 = ax.bar(x, esm, w, label="ESM-MLP (sequence)", color="#2ca02c",
                edgecolor="white", zorder=3)
    b3 = ax.bar(x + w, clf, w, label="ReliFuse (fusion)", color="#d62728",
                edgecolor="white", zorder=3)
    for bars in (b1, b2, b3):
        for r in bars:
            ax.text(r.get_x() + r.get_width() / 2, r.get_height() + 0.008,
                    f"{r.get_height():.3f}", ha="center", va="bottom",
                    fontsize=10, fontweight="bold")
    for i, d in enumerate(deltas):
        ax.text(x[i] + w, clf[i] - 0.028, d, ha="center", va="top",
                fontsize=9.5, color="#d62728", fontstyle="italic",
                bbox=dict(boxstyle="round,pad=0.25", fc="#fff0f0",
                          ec="#d62728", alpha=0.85))
    ax.set_xticks(x)
    ax.set_xticklabels(strata)
    ax.set_ylabel("MCC")
    ax.set_ylim(0.5, 0.88)
    ax.set_title("Confidence-conditional complementarity:\n"
                 "three-way comparison by pLDDT stratum (20-seed aggregate)")
    ax.legend(loc="upper right", ncol=1, framealpha=0.95)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig4_plddt_threeway.png"), dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved fig4_plddt_threeway.png")


# ---------------------------------------------------------------- bias b
def fig_bias():
    df = load_agg()
    from scipy import stats as sps
    b = df["b_gated2"].to_numpy(float)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.hist(b, bins=70, density=True, color="#9467bd", alpha=0.55,
            label="density histogram")
    xs = np.linspace(b.min(), b.max(), 300)
    ax.plot(xs, sps.gaussian_kde(b)(xs), color="#4b2e83", lw=2, label="KDE")
    for v, c, lab in [(b.mean(), "#d62728", "mean"),
                      (np.median(b), "#1f77b4", "median")]:
        ax.axvline(v, color=c, ls="--", lw=1.2)
        ax.text(v, ax.get_ylim()[1] * 0.94, f"{lab} {v:.3f}", rotation=90,
                va="top", ha="right", fontsize=8, color=c)
    ax.set_xlabel("learned logit bias b")
    ax.set_ylabel("density")
    ax.set_title("ReliFuse bias distribution")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)

    p1 = df["p_mkmp"].to_numpy(float)
    p2 = df["p_gated2"].to_numpy(float)
    disagree = np.abs(p1 - p2)
    ax2 = axes[1]
    tert = np.quantile(disagree, [0, 1 / 3, 2 / 3, 1])
    labels = ["low disagreement", "mid disagreement", "high disagreement"]
    groups = [b[(disagree >= lo) & (disagree < hi)]
              for lo, hi in zip(tert[:-1], tert[1:])]
    groups[-1] = b[disagree >= tert[-2]]
    parts = ax2.violinplot(groups, showmeans=True, showmedians=True)
    for pc in parts["bodies"]:
        pc.set_facecolor("#9467bd")
        pc.set_alpha(0.65)
    ax2.set_xticks([1, 2, 3])
    ax2.set_xticklabels(labels, rotation=12, ha="right")
    ax2.set_ylabel("learned logit bias b")
    ax2.set_title("bias b by model disagreement")
    ax2.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_bias_distribution.png"), dpi=300)
    plt.close(fig)
    print("saved fig_bias_distribution.png")


# ---------------------------------------------------------------- 0924 figs
C_MKMP, C_LGF = "#4c72b0", "#55a868"
METHOD_CN = {"mkmp": "MKMP", "weighted": "weighted", "lgf": "ReliFuse"}


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_e2():
    df = pd.read_csv(os.path.join(RES0924, "e2_overlap_strata_metrics.csv"))
    df = df[df.protocol == "finetuned"]
    subsets = ["full", "accession_disjoint", "seq30_disjoint",
               "seq50_disjoint", "no_balf_overlap"]
    labels = {"full": "full", "accession_disjoint": "accession-disjoint",
              "seq30_disjoint": "seq30-disjoint",
              "seq50_disjoint": "seq50-disjoint",
              "no_balf_overlap": "no-BALF-overlap"}
    x = np.arange(len(subsets))
    w = 0.36
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    for j, m in enumerate(["mkmp", "lgf"]):
        vals = [df[(df.subset == s) & (df.method == m)]["MCC"].iloc[0]
                for s in subsets]
        ax.bar(x + (j - 0.5) * w, vals, w, label=METHOD_CN[m],
               color=C_MKMP if m == "mkmp" else C_LGF)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[s] for s in subsets], fontsize=8)
    ax.set_ylabel("MCC")
    ax.set_title("E2  GlycoHBF 重叠/同源分层（finetuned）")
    ax.legend(frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_e2_overlap_strata.png"),
                bbox_inches="tight")
    plt.close(fig)
    print("saved fig_e2_overlap_strata.png")


def fig_e8():
    df = pd.read_csv(os.path.join(RES0924, "e8_label_sensitivity_metrics.csv"))
    order = ["full", "no_plasma_only", "no_intracellular",
             "extracellular_only", "no_contaminants", "secreted_strict"]
    labels = {"full": "full", "no_plasma_only": "no-plasma",
              "no_intracellular": "no-intracellular",
              "extracellular_only": "extracellular",
              "no_contaminants": "no-contaminant",
              "secreted_strict": "secreted-strict"}
    x = np.arange(len(order))
    w = 0.36
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    for j, m in enumerate(["mkmp", "lgf"]):
        vals = [df[(df.stratum == s) & (df.method == m)]["MCC"].iloc[0]
                for s in order]
        ax.bar(x + (j - 0.5) * w, vals, w, label=METHOD_CN[m],
               color=C_MKMP if m == "mkmp" else C_LGF)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[s] for s in order], fontsize=8)
    ax.set_ylabel("MCC")
    ax.set_title("E8  严格标签敏感性")
    ax.legend(frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_e8_label_sensitivity.png"),
                bbox_inches="tight")
    plt.close(fig)
    print("saved fig_e8_label_sensitivity.png")


# ---------------------------------------------------------------- GO
def _load_obo(path):
    term, cur = {}, None
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                cur = {}
            elif line == "" and cur is not None:
                if "id" in cur and "name" in cur:
                    term[cur["id"]] = (cur["name"], cur.get("namespace", ""))
                cur = None
            elif cur is not None and ": " in line:
                k, v = line.split(": ", 1)
                cur[k] = v
    return term


def _load_gaf(path):
    acc = defaultdict(set)
    import gzip
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 9 and "NOT" not in p[3]:
                acc[p[1]].add(p[4])
    return acc


def _top500(prob_col):
    files = sorted(glob.glob(os.path.join(
        V1, "results", "logit_fusion", "preds_esm150m_v4_s*.csv")))
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    g = df.groupby("name").agg(y=("y", "max"), p=(prob_col, "mean")).reset_index()
    return set(g.sort_values("p", ascending=False)["name"].head(500))


def _enrich(fg, bg, obo):
    fg = {a for a in fg if a in bg}
    all_go = sorted({g for a in bg for g in bg[a]})
    N = len(bg)
    K_g = defaultdict(int)
    for a in bg:
        for g in bg[a]:
            K_g[g] += 1
    rows = []
    n = len(fg)
    for go in all_go:
        k = K_g[go]
        x = sum(1 for a in fg if go in bg[a])
        if k == 0 or x == 0:
            continue
        a_, b_ = x, k - x
        c_, d_ = n - x, (N - k) - (n - x)
        if min(a_, b_, c_, d_) < 0:
            continue
        _, p = fisher_exact([[a_, b_], [c_, d_]], alternative="greater")
        rows.append((go, k, x, n, p))
    df = pd.DataFrame(rows, columns=["GO", "bg", "fg", "n_fg", "p"])
    df = df[df["fg"] >= 2].copy()
    df["FDR"] = (df["p"].values * len(df) / df["p"].rank(method="first").values)
    df["FDR"] = df["FDR"].clip(upper=1)
    df = df[df["FDR"] < 0.05]
    df["name"] = df["GO"].map(lambda g: obo.get(g, (g, ""))[0])
    return df.sort_values("FDR")


def fig_go():
    obo = _load_obo(os.path.join(ROOT, "go_enrichment", "go-basic.obo"))
    bg = _load_gaf(os.path.join(ROOT, "go_enrichment", "goa_human.gaf"))
    lgf = _top500("p_gated2")
    res = _enrich(lgf, bg, obo)
    res.to_csv(os.path.join(V1, "results", "go_enrichment",
                            "go_enrich_lgf.csv"), index=False)
    r1 = res.head(12)
    labels = [f"{n} ({g})" for n, g in zip(r1["name"], r1["GO"])]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(np.arange(len(labels))[::-1],
            -np.log10(r1["FDR"].clip(lower=1e-12)), color="#2ca02c")
    ax.set_yticks(np.arange(len(labels))[::-1])
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("-log10 FDR")
    ax.set_title("ReliFuse Top-500 GO enrichment")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_go_enrichment.png"), dpi=300)
    plt.close(fig)
    print("saved fig_go_enrichment.png")


def fig_safety_switch():
    """Re-render with the paper's v44/v45 terminology (consistency gating)."""
    df = pd.read_csv(os.path.join(RES0924, "structure_safety_switch.csv"))
    metrics = ["seq_auc", "af_auc", "ef_auc", "gated_auc"]
    labels = ["ESM-MLP\n(sequence)", "AlphaFold\nstructure",
              "ESMFold\nstructure", "Consistency\ngate"]
    means = [df[m].mean() for m in metrics]
    sds = [df[m].std(ddof=1) for m in metrics]
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    x = np.arange(len(metrics))
    ax.bar(x, means, yerr=sds, capsize=5,
           color=["#4c72b0", "#dd8452", "#c44e52", "#55a868"])
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("AUC (mean +/- SD)")
    ax.set_ylim(0, 1.0)
    ax.set_title("Structure-source consistency gating (n=86)")
    for i, v in enumerate(means):
        ax.text(i, v + sds[i] + 0.02, f"{v:.3f}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_safety_switch.png"), dpi=200)
    plt.close(fig)
    print("saved fig_safety_switch.png")


if __name__ == "__main__":
    fig4()
    fig_bias()
    fig_e2()
    fig_e8()
    fig_go()
    fig_safety_switch()
