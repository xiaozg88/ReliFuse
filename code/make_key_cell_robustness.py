# -*- coding: utf-8 -*-
"""Statistical robustness of the double-uncertainty key-cell fusion gain.

Key cell: pLDDT < 70 AND maximum test-to-train sequence identity < 0.20
(n = 273, 22 positives).  The point estimate of the fusion gain is
dMCC = MCC(ReliFuse) - max(MCC(MKMP), MCC(ESM-MLP)) = +0.1755 on the
20-seed pooled benchmark predictions.  Because the cell contains only 22
positives, this script adds three checks requested for reviewer-proofing:

1. Cluster/bootstrap 95% CI (10,000 resamples, seed 3407).
2. Label-permutation two-sided p-value (10,000 permutations, seed 3407):
   shuffle the detection labels within the cell; this preserves the cell's
   class composition while breaking the label-prediction association, so it
   tests whether the observed dMCC exceeds chance given the cell alone.
3. Leave-one-positive-out (LOPO) sensitivity: drop each of the 22 positives
   in turn and recompute dMCC, checking whether a single positive drives
   the gain.

Output: results/e25_key_cell_robustness.csv
"""

import os

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
DATADIR = os.path.join(ROOT, "data")
RESDIR = os.path.join(ROOT, "results")

N_BOOT = 10000
N_PERM = 10000
SEED = 3407


def load_key_cell():
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
    cell = df[(df["plddt"] < 70.0) & (df["max_identity"] < 0.20)].reset_index(drop=True)
    return cell


def dmcc(y, p_mkmp, p_esm, p_fuse):
    m = matthews_corrcoef(y, (p_mkmp >= 0.5).astype(int))
    e = matthews_corrcoef(y, (p_esm >= 0.5).astype(int))
    f = matthews_corrcoef(y, (p_fuse >= 0.5).astype(int))
    return f - max(m, e)


def main():
    cell = load_key_frame = load_key_cell()
    y = cell["y"].values.astype(int)
    p_mkmp = cell["p_mkmp"].values
    p_esm = cell["p_esm"].values
    p_fuse = cell["p_fuse"].values
    n, n_pos = len(cell), int(y.sum())
    assert (n, n_pos) == (273, 22), f"unexpected cell size {(n, n_pos)}"

    obs = dmcc(y, p_mkmp, p_esm, p_fuse)

    rng = np.random.default_rng(SEED)

    # 1) bootstrap 95% CI
    boot = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = rng.choice(n, size=n, replace=True)
        yy = y[idx]
        if yy.sum() in (0, len(yy)):          # degenerate draw
            boot[b] = np.nan
            continue
        boot[b] = dmcc(yy, p_mkmp[idx], p_esm[idx], p_fuse[idx])
    boot = boot[~np.isnan(boot)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    frac_pos = float((boot > 0).mean())

    # 2) within-cell label permutation (two-sided)
    perm = np.empty(N_PERM)
    for b in range(N_PERM):
        yy = rng.permutation(y)
        perm[b] = dmcc(yy, p_mkmp, p_esm, p_fuse)
    p_two = float((np.abs(perm) >= abs(obs)).mean())

    # 3) leave-one-positive-out
    pos_idx = np.flatnonzero(y == 1)
    lopo = np.empty(len(pos_idx))
    for k, j in enumerate(pos_idx):
        keep = np.ones(n, dtype=bool)
        keep[j] = False
        lopo[k] = dmcc(y[keep], p_mkmp[keep], p_esm[keep], p_fuse[keep])

    row = {
        "cell": "pLDDT<70 & identity<0.20",
        "n": n,
        "n_pos": n_pos,
        "point_dMCC": round(float(obs), 4),
        "boot_n_resamples": int(len(boot)),
        "boot_CI95_low": round(float(lo), 4),
        "boot_CI95_high": round(float(hi), 4),
        "boot_frac_positive": round(frac_pos, 4),
        "perm_n": N_PERM,
        "perm_p_two_sided": round(p_two, 5),
        "lopo_folds": len(pos_idx),
        "lopo_min": round(float(lopo.min()), 4),
        "lopo_max": round(float(lopo.max()), 4),
        "lopo_mean": round(float(lopo.mean()), 4),
        "lopo_all_positive": bool((lopo > 0).all()),
    }
    out = pd.DataFrame([row])
    os.makedirs(RESDIR, exist_ok=True)
    path = os.path.join(RESDIR, "e25_key_cell_robustness.csv")
    out.to_csv(path, index=False)
    print(f"saved {path}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
