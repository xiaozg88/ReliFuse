# -*- coding: utf-8 -*-
"""Reliability-gated complementarity: table data + continuous characterization.

Outputs
-------
- ``results/e20_two_axis_complementarity.csv``
  The 3x3 (pLDDT x sequence-identity) grid of MCC / AUC / dMCC / dAUC, which
  feeds Table I and Table II of the paper.
- ``figures/fig_two_axis_complementarity.png/.pdf``
  A continuous, non-redundant companion to the discrete tables:
    (a) MCC vs pLDDT for MKMP / ESM-MLP / ReliFuse (the confidence crossover);
    (b) fusion gain dMCC vs pLDDT, split by homology-free (<0.20) vs the rest.
"""

import os

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
DATADIR = os.path.join(ROOT, "data")
RESDIR = os.path.join(ROOT, "results")
FIGDIR = os.path.join(ROOT, "figures")

# Okabe-Ito colour-blind-safe palette
C_MKMP = "#0072B2"   # blue  (structure pathway)
C_ESM = "#E69F00"    # orange (sequence pathway)
C_FUSE = "#009E73"   # green  (fusion)
C_LOID = "#D55E00"   # vermilion (homology-free)
C_HIID = "#56B4E9"   # sky blue (homology-rich)


def load_test_frame():
    mkmp = pd.read_csv(os.path.join(DATADIR, "benchmark_pred_MKMP.csv"))[
        ["name", "y", "p"]].rename(columns={"p": "p_mkmp"})
    esm = pd.read_csv(os.path.join(DATADIR, "benchmark_pred_ESM_MLP.csv"))[
        ["name", "p"]].rename(columns={"p": "p_esm"})
    fuse = pd.read_csv(os.path.join(DATADIR, "benchmark_pred_ReliFuse.csv"))[
        ["name", "p"]].rename(columns={"p": "p_fuse"})
    hom = pd.read_csv(os.path.join(DATADIR, "homology_identity.csv"))[
        ["name", "max_identity"]]

    pn = np.load(os.path.join(DATADIR, "plddt_names_8069.npy"), allow_pickle=True)
    pf = np.load(os.path.join(DATADIR, "plddt_features_8069.npy"), allow_pickle=True)
    pmap = {str(n): float(pf[i, 0]) for i, n in enumerate(pn)}

    df = mkmp.merge(esm, on="name").merge(fuse, on="name").merge(hom, on="name")
    df["plddt"] = df["name"].map(pmap)
    df = df.dropna(subset=["plddt"]).reset_index(drop=True)
    return df


def _mcc(y, p):
    return matthews_corrcoef(y, (p >= 0.5).astype(int))


def cell_metrics(g):
    y = g["y"].values
    out = {"n": len(g), "n_pos": int(y.sum()), "pos_rate": float(y.mean())}
    for col in ("p_mkmp", "p_esm", "p_fuse"):
        out[col + "_mcc"] = _mcc(y, g[col].values)
        out[col + "_auc"] = roc_auc_score(y, g[col].values)
    best = "MKMP" if out["p_mkmp_mcc"] >= out["p_esm_mcc"] else "ESM"
    out["best_single"] = best
    out["best_mcc"] = max(out["p_mkmp_mcc"], out["p_esm_mcc"])
    out["dMCC"] = out["p_fuse_mcc"] - out["best_mcc"]
    out["dAUC"] = out["p_fuse_auc"] - max(out["p_mkmp_auc"], out["p_esm_auc"])
    return out


def make_grid_csv(df):
    plddt_bins = [(-1.0, 70.0, "<70"), (70.0, 90.0, "70-90"), (90.0, 101.0, ">90")]
    id_bins = [(-0.001, 0.20, "<0.20"), (0.20, 0.40, "0.20-0.40"), (0.40, 1.01, ">0.40")]
    rows = []
    for plo, phi, plab in plddt_bins:
        for ilo, ihi, ilab in id_bins:
            m = (df["plddt"] >= plo) & (df["plddt"] < phi)
            m &= (df["max_identity"] >= ilo) & (df["max_identity"] < ihi)
            g = df[m]
            if len(g) == 0:
                continue
            r = cell_metrics(g)
            rows.append({"plddt_stratum": plab, "identity_stratum": ilab,
                         **{k: (round(v, 4) if isinstance(v, (int, float)) else v)
                            for k, v in r.items()}})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(RESDIR, "e20_two_axis_complementarity.csv"), index=False)
    return res


def axis_orthogonality(df):
    """Spearman/Pearson correlation between the two reliability axes on the
    test set (quoted in Section IV-A of the main text).  Added in the v51
    traceability pass so the in-text statistics are produced by a shipped
    script rather than computed offline."""
    from scipy import stats as sps
    sp = sps.spearmanr(df["plddt"], df["max_identity"])
    pe = sps.pearsonr(df["plddt"], df["max_identity"])
    out = pd.DataFrame([{
        "n": len(df),
        "spearman_r": round(float(sp.statistic), 4),
        "spearman_p": float(sp.pvalue),
        "pearson_r": round(float(pe.statistic), 4),
        "pearson_p": float(pe.pvalue),
    }])
    out.to_csv(os.path.join(RESDIR, "e20_axis_orthogonality.csv"), index=False)
    print("axis orthogonality: Spearman r = %.4f, Pearson r = %.4f (n=%d)"
          % (sp.statistic, pe.statistic, len(df)))
    return out


def make_continuous_figure(df):
    # equal-frequency pLDDT bins (continuous x-axis)
    df = df.sort_values("plddt").reset_index(drop=True)
    df["bin"] = pd.qcut(df["plddt"], 8, labels=False, duplicates="drop")
    nbins = int(df["bin"].max()) + 1

    # panel (a): MCC vs pLDDT for the three models
    a_x, a_mkmp, a_esm, a_fuse = [], [], [], []
    # panel (b): dMCC vs pLDDT split by homology-free vs rest
    b_x_lo, b_lo = [], []
    b_x_hi, b_hi = [], []
    for b in range(nbins):
        g = df[df["bin"] == b]
        if len(g) == 0:
            continue
        x = float(g["plddt"].mean())
        a_x.append(x)
        a_mkmp.append(_mcc(g["y"], g["p_mkmp"]))
        a_esm.append(_mcc(g["y"], g["p_esm"]))
        a_fuse.append(_mcc(g["y"], g["p_fuse"]))
        # dMCC within homology-free (<0.20) and homology-rich (>=0.20)
        for grp, xarr, yarr in ((g[g["max_identity"] < 0.20], b_x_lo, b_lo),
                                (g[g["max_identity"] >= 0.20], b_x_hi, b_hi)):
            if len(grp) == 0:
                continue
            best = max(_mcc(grp["y"], grp["p_mkmp"]), _mcc(grp["y"], grp["p_esm"]))
            dm = _mcc(grp["y"], grp["p_fuse"]) - best
            xarr.append(x)
            yarr.append(dm)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1))

    ax = axes[0]
    ax.plot(a_x, a_mkmp, "-o", ms=4, lw=1.4, color=C_MKMP, label="MKMP (structure)")
    ax.plot(a_x, a_esm, "-s", ms=4, lw=1.4, color=C_ESM, label="ESM-MLP (sequence)")
    ax.plot(a_x, a_fuse, "-^", ms=4, lw=1.4, color=C_FUSE, label="ReliFuse (fusion)")
    ax.set_xlabel("Mean pLDDT (bin)")
    ax.set_ylabel("MCC")
    ax.set_title("(a) Three-way MCC vs pLDDT", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(alpha=0.25)

    ax = axes[1]
    ax.axhline(0, color="0.6", lw=0.8, ls="--")
    ax.plot(b_x_lo, b_lo, "-o", ms=4, lw=1.4, color=C_LOID,
            label="identity < 0.20")
    ax.plot(b_x_hi, b_hi, "-s", ms=4, lw=1.4, color=C_HIID,
            label="identity \u2265 0.20")
    ax.set_xlabel("Mean pLDDT (bin)")
    ax.set_ylabel("$\\Delta$MCC (ReliFuse $-$ best single)")
    ax.set_title("(b) Fusion gain vs pLDDT, by homology", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_two_axis_complementarity.png"), dpi=400)
    fig.savefig(os.path.join(FIGDIR, "fig_two_axis_complementarity.pdf"))
    plt.close(fig)


def main():
    df = load_test_frame()
    assert len(df) == 1613, f"expected 1613 test proteins, got {len(df)}"
    make_grid_csv(df)
    axis_orthogonality(df)
    make_continuous_figure(df)
    print("saved", os.path.join(RESDIR, "e20_two_axis_complementarity.csv"))
    print("saved figure", os.path.join(FIGDIR, "fig_two_axis_complementarity.png"))


if __name__ == "__main__":
    main()
