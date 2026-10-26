# -*- coding: utf-8 -*-
"""Two-axis robustness (2x2) and explicit reliability-gate ablation."""
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from sklearn.metrics import matthews_corrcoef, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
DATADIR = os.path.join(ROOT, "data")
RESDIR = os.path.join(ROOT, "results")
SEEDS = [3407, 42, 123, 2024, 777, 1, 2, 3, 4, 5, 6, 7, 8, 10, 21, 99, 202, 555, 888, 1000]


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -50, 50)))


def mcc(y, p):
    return matthews_corrcoef(y, (p >= 0.5).astype(int))


def load_probability_pooled():
    mkmp = pd.read_csv(os.path.join(DATADIR, "benchmark_pred_MKMP.csv"))[["name", "y", "p"]].rename(columns={"p": "p_mkmp"})
    esm = pd.read_csv(os.path.join(DATADIR, "benchmark_pred_ESM_MLP.csv"))[["name", "p"]].rename(columns={"p": "p_esm"})
    fuse = pd.read_csv(os.path.join(DATADIR, "benchmark_pred_ReliFuse.csv"))[["name", "p"]].rename(columns={"p": "p_fuse"})
    hom = pd.read_csv(os.path.join(DATADIR, "homology_identity.csv"))[["name", "max_identity"]]
    pn = np.load(os.path.join(DATADIR, "plddt_names_8069.npy"), allow_pickle=True)
    pf = np.load(os.path.join(DATADIR, "plddt_features_8069.npy"), allow_pickle=True)
    pmap = {str(n): float(pf[i, 0]) for i, n in enumerate(pn)}
    df = mkmp.merge(esm, on="name").merge(fuse, on="name").merge(hom, on="name")
    df["plddt"] = df["name"].map(pmap)
    return df.dropna(subset=["plddt"]).reset_index(drop=True)


def load_logit_pooled():
    frames = [pd.read_csv(os.path.join(DATADIR, f"preds_esm150m_v4_s{s}.csv"))[["name", "z_mkmp_cal", "z_esm_cal", "g_gated2", "b_gated2"]] for s in SEEDS]
    agg = pd.concat(frames).groupby("name").mean().reset_index()
    hom = pd.read_csv(os.path.join(DATADIR, "homology_identity.csv"))[["name", "y", "max_identity"]]
    pn = np.load(os.path.join(DATADIR, "plddt_names_8069.npy"), allow_pickle=True)
    pf = np.load(os.path.join(DATADIR, "plddt_features_8069.npy"), allow_pickle=True)
    pmap = {str(n): float(pf[i, 0]) for i, n in enumerate(pn)}
    df = agg.merge(hom, on="name")
    df["plddt"] = df["name"].map(pmap)
    df = df.dropna(subset=["plddt"]).reset_index(drop=True)
    z = df.z_mkmp_cal + df.g_gated2 * (df.z_esm_cal - df.z_mkmp_cal) + df.b_gated2
    df["p_learned"] = sigmoid(z)
    return df


def cell_pp(df, plo, phi, ilo, ihi):
    g = df[(df.plddt >= plo) & (df.plddt < phi) & (df.max_identity >= ilo) & (df.max_identity < ihi)]
    best_mcc = max(mcc(g.y, g.p_mkmp), mcc(g.y, g.p_esm))
    best_auc = max(roc_auc_score(g.y, g.p_mkmp), roc_auc_score(g.y, g.p_esm))
    return {
        "n": len(g), "n_pos": int(g.y.sum()),
        "dMCC": mcc(g.y, g.p_fuse) - best_mcc,
        "dAUC": roc_auc_score(g.y, g.p_fuse) - best_auc,
    }


def main():
    pp = load_probability_pooled()

    # #3: 2x2 robustness (probability-pooled)
    rows = []
    for plo, phi, plab in [(-1, 70, "<70"), (70, 101, ">=70")]:
        for ilo, ihi, ilab in [(-0.001, 0.20, "<0.20"), (0.20, 1.01, ">=0.20")]:
            r = cell_pp(pp, plo, phi, ilo, ihi)
            rows.append({"plddt": plab, "identity": ilab,
                         **{k: (round(v, 4) if isinstance(v, (int, float)) else v) for k, v in r.items()}})
    pd.DataFrame(rows).to_csv(os.path.join(RESDIR, "e21_two_by_two_robustness.csv"), index=False)
    print("== 2x2 robustness (probability-pooled) ==")
    print(pd.DataFrame(rows).pivot(index="plddt", columns="identity", values="dMCC").round(4).to_string())

    # #4: explicit reliability gate + learned-gate correlation
    df = load_logit_pooled()
    corr_rows = []
    for col in ["g_gated2", "b_gated2"]:
        for feat in ["plddt", "max_identity"]:
            r, p = stats.spearmanr(df[col], df[feat])
            corr_rows.append({"component": col, "feature": feat, "spearman_r": round(r, 4), "p": f"{p:.2e}"})

    z1, z2, y = df.z_mkmp_cal.values, df.z_esm_cal.values, df.y.values
    xp = (df.plddt.values - df.plddt.mean()) / df.plddt.std()
    xi = (df.max_identity.values - df.max_identity.mean()) / df.max_identity.std()
    X = np.stack([xp, xi, np.ones_like(xp)], 1)

    def loss(theta):
        g = sigmoid(X @ theta[:3])
        b = X @ theta[3:]
        p = sigmoid(z1 + g * (z2 - z1) + b)
        return -np.mean(y * np.log(p + 1e-12) + (1 - y) * np.log(1 - p + 1e-12))

    res = minimize(loss, np.zeros(6), method="L-BFGS-B")
    g_ex, b_ex = sigmoid(X @ res.x[:3]), X @ res.x[3:]
    df["p_explicit"] = sigmoid(z1 + g_ex * (z2 - z1) + b_ex)

    def key_cell(pcol):
        m = df[(df.plddt < 70) & (df.max_identity < 0.20)]
        bm = max(mcc(m.y, sigmoid(m.z_mkmp_cal)), mcc(m.y, sigmoid(m.z_esm_cal)))
        return mcc(m.y, m[pcol]) - bm

    summary = {
        "g_explicit_mean": round(float(g_ex.mean()), 4), "g_explicit_sd": round(float(g_ex.std()), 4),
        "b_explicit_mean": round(float(b_ex.mean()), 4), "b_explicit_sd": round(float(b_ex.std()), 4),
        "overall_mcc_learned": round(mcc(df.y, df.p_learned), 4), "overall_mcc_explicit": round(mcc(df.y, df.p_explicit), 4),
        "overall_auc_learned": round(roc_auc_score(df.y, df.p_learned), 4), "overall_auc_explicit": round(roc_auc_score(df.y, df.p_explicit), 4),
        "key_cell_dMCC_learned": round(key_cell("p_learned"), 4),
        "key_cell_dMCC_explicit": round(key_cell("p_explicit"), 4),
    }

    pd.DataFrame(corr_rows).to_csv(os.path.join(RESDIR, "e22_gate_correlations.csv"), index=False)
    pd.DataFrame([summary]).to_csv(os.path.join(RESDIR, "e23_explicit_gate_summary.csv"), index=False)

    print("\n== learned-gate correlations (Spearman) ==")
    print(pd.DataFrame(corr_rows).to_string(index=False))
    print("\n== explicit-gate summary ==")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
