# -*- coding: utf-8 -*-
"""Validation-fitted explicit reliability gate (strict upgrade of #4b).

The learned gate in ReliFuse is a small MLP over ~16 features.  This script
replaces it with an *explicit* affine gate that uses ONLY the two reliability
axes (structural confidence pLDDT and evolutionary support sequence identity),
fits its six parameters on the fixed VALIDATION split, and evaluates once on
the TEST split -- the methodologically correct ablation that removes the
in-test-set-fitting optimism of the earlier diagnostic.

Run from the repo root (or anywhere)::

    python code/make_explicit_gate_valid.py --seeds 3407,42,123,2024,777

It reuses the base-model logits already computed by
``v1-01/v1/run_logit_fusion_pilot.py`` (MKMP v12 checkpoints + the ESM mean-pool
cache) and, for each seed, (1) obtains calibrated structure/sequence logits on
validation and test, (2) caches them to CSV, then (3) fits the explicit gate on
validation and reports overall MCC/AUC and the two-axis dMCC map on test.

Note on the identity feature: the exact test-to-train maximum sequence identity
is only precomputed for the TEST set (``homology_identity_test.csv``).  For the
VALIDATION set we use a coarse cluster-overlap proxy from
``cluster_folds_50pct.csv`` (same 50% CD-HIT cluster as a training protein ->
0.5, else 0.0); replace ``load_homology`` with MMseqs2 val-to-train identity
for an exact version.
"""

import argparse
import csv
import os
import sys

import numpy as np
import torch
from scipy.optimize import minimize
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from torch_geometric.loader import DataLoader

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
RESDIR = os.path.join(ROOT, "results")
os.makedirs(RESDIR, exist_ok=True)

# Reuse the v1-01 pipeline. This script depends on the upstream training
# pipeline (run_logit_fusion_pilot, pLDDT features), which lives four levels
# above this file in the development tree and is NOT part of the release
# repository; the released per-protein predictions under ../data/ are the
# supported entry point for reproduction.
V1 = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "v1"))
if V1 not in sys.path:
    sys.path.insert(0, V1)
sys.path.insert(0, os.path.normpath(os.path.join(V1, "..", "Mgate")))

import run_logit_fusion_pilot as rlfp  # noqa: E402

SPLIT = os.path.join(ROOT, "data", "splits", "main_split_seed3407.csv")
TEST_HOM = os.path.join(ROOT, "data", "splits", "homology_identity_test.csv")
CLUSTER = os.path.join(ROOT, "data", "splits", "cluster_folds_50pct.csv")
PLDDT_FEAT = os.path.join(V1, "data", "plddt_features_8069.npy")
PLDDT_NAMES = os.path.join(V1, "data", "plddt_names_8069.npy")


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(x, float), -50, 50)))


def load_plddt_map():
    feats = np.load(PLDDT_FEAT)
    names = np.load(PLDDT_NAMES, allow_pickle=True)
    return {str(n): float(feats[i, 0]) for i, n in enumerate(names)}


def load_homology(val_names):
    """Return name -> evolutionary-support score (test exact, val proxy)."""
    hom = {}
    import pandas as pd

    # test: exact max test-to-train identity
    if os.path.exists(TEST_HOM):
        df = pd.read_csv(TEST_HOM)
        hom.update(dict(zip(df["name"].astype(str), df["max_identity"].astype(float))))

    # validation: cluster-overlap proxy (0.5 if shared 50% cluster w/ training)
    if os.path.exists(SPLIT) and os.path.exists(CLUSTER):
        split = pd.read_csv(SPLIT)
        train = set(split.loc[split["split"] == "train", "name"].astype(str))
        cl = pd.read_csv(CLUSTER)
        train_clusters = set(cl.loc[cl["name"].astype(str).isin(train), "cluster_id"])
        name_to_cluster = dict(zip(cl["name"].astype(str), cl["cluster_id"]))
        for n in val_names:
            if n not in hom:
                cid = name_to_cluster.get(n)
                hom[n] = 0.5 if cid in train_clusters else 0.0
    return hom


def explicit_gate_fit(z1, z2, y, plddt, ident):
    xp = (np.asarray(plddt) - np.mean(plddt)) / (np.std(plddt) + 1e-8)
    xi = (np.asarray(ident) - np.mean(ident)) / (np.std(ident) + 1e-8)
    X = np.stack([xp, xi, np.ones_like(xp)], axis=1)
    z1, z2, y = np.asarray(z1), np.asarray(z2), np.asarray(y)

    def loss(theta):
        g = sigmoid(X @ theta[:3])
        b = X @ theta[3:]
        p = sigmoid(z1 + g * (z2 - z1) + b)
        p = np.clip(p, 1e-7, 1 - 1e-7)
        return -float(np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))

    res = minimize(loss, np.zeros(6), method="L-BFGS-B")
    return res.x


def apply_explicit(theta, z1, z2, plddt, ident, p_mean, i_mean, p_std, i_std):
    xp = (np.asarray(plddt) - p_mean) / (p_std + 1e-8)
    xi = (np.asarray(ident) - i_mean) / (i_std + 1e-8)
    X = np.stack([xp, xi, np.ones_like(xp)], axis=1)
    g = sigmoid(X @ theta[:3])
    b = X @ theta[3:]
    return sigmoid(np.asarray(z1) + g * (np.asarray(z2) - np.asarray(z1)) + b)


def mcc(y, p):
    return matthews_corrcoef(y, (p >= 0.5).astype(int))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="3407,42,123,2024,777")
    ap.add_argument("--epochs", type=int, default=500)
    ap.add_argument("--patience", type=int, default=20)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--esm-cache", default=rlfp.ESM_CACHE)
    ap.add_argument("--val-cache", default=os.path.join(RESDIR, "val_logits_cache.csv"))
    ap.add_argument("--skip-inference", action="store_true")
    args = ap.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    device = torch.device(args.device)
    plddt_map = load_plddt_map()

    if not args.skip_inference:
        data_list, idx, dataset = rlfp.load_mkmp_dataset()
        val_loader = DataLoader(
            torch.utils.data.Subset(dataset, idx["val"]), batch_size=64, shuffle=False)
        test_loader = DataLoader(
            torch.utils.data.Subset(dataset, idx["test"]), batch_size=64, shuffle=False)
        esm_names, esm_X, esm_y, esm_splits = rlfp.load_esm_cache(args.esm_cache)

        cache_rows = []
        for seed in seeds:
            mk_path = rlfp.MKMP_MODEL.format(seed=seed)
            rlfp.reseed(seed)
            mk_model = rlfp.build_mkmp_model(dataset, device)
            mk_model.load_state_dict(torch.load(mk_path, weights_only=True, map_location="cpu"))
            mk_model.to(device)
            mk_val_z, mk_val_y = rlfp.eval_mkmp(mk_model, val_loader, device)
            mk_test_z, mk_test_y = rlfp.eval_mkmp(mk_model, test_loader, device)
            del mk_model

            esm_val_z, esm_val_y, esm_test_z, esm_test_y = rlfp.train_esm(
                seed, esm_names, esm_X, esm_y, esm_splits, device,
                args.epochs, args.patience, rlfp.Logger(os.devnull))

            val_names, val_y, val_z1, val_z2 = rlfp.align_dicts(
                mk_val_z, mk_val_y, esm_val_z, esm_val_y)
            test_names, test_y, test_z1, test_z2 = rlfp.align_dicts(
                mk_test_z, mk_test_y, esm_test_z, esm_test_y)

            t1 = rlfp.fit_temperature(val_z1, val_y)
            t2 = rlfp.fit_temperature(val_z2, val_y)
            for split, names, z1, z2 in [("val", val_names, val_z1 / t1, val_z2 / t2),
                                         ("test", test_names, test_z1 / t1, test_z2 / t2)]:
                for i, n in enumerate(names):
                    cache_rows.append([seed, split, n, z1[i], z2[i], plddt_map.get(n, 70.0)])

        with open(args.val_cache, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["seed", "split", "name", "z_mkmp_cal", "z_esm_cal", "plddt"])
            w.writerows(cache_rows)
        print("saved logits cache:", args.val_cache, "rows=", len(cache_rows))

    import pandas as pd
    cache = pd.read_csv(args.val_cache)
    cache["name"] = cache["name"].astype(str)
    val_names = cache.loc[cache["split"] == "val", "name"].unique().tolist()
    hom = load_homology(val_names)

    # Fit the explicit gate on the (freshly cached) validation logits only.
    val = cache[cache["split"] == "val"].groupby("name").agg(
        z_mkmp_cal=("z_mkmp_cal", "mean"), z_esm_cal=("z_esm_cal", "mean"),
        plddt=("plddt", "mean")).reset_index()
    val["identity"] = val["name"].map(hom)
    # labels for val
    import run_logit_fusion_pilot as _r
    data_list, idx, _ = _r.load_mkmp_dataset()
    val_name_to_y = {data_list[i].name: int(data_list[i].y) for i in idx["val"]}
    val["y"] = val["name"].map(val_name_to_y)
    val = val.dropna(subset=["y", "identity"])
    theta = explicit_gate_fit(val.z_mkmp_cal, val.z_esm_cal, val.y, val.plddt, val.identity)
    print("explicit gate theta (g: a0,a1,a2 ; b: b0,b1,b2) =", np.round(theta, 4))

    # Evaluate both gates on the SAME official test logits (preds files).
    lf = []
    for seed in seeds:
        df = pd.read_csv(os.path.join(ROOT, "data", f"preds_esm150m_v4_s{seed}.csv"))
        lf.append(df[["name", "z_mkmp_cal", "z_esm_cal", "g_gated2", "b_gated2"]])
    lf = pd.concat(lf).groupby("name").mean().reset_index()
    z = lf.z_mkmp_cal + lf.g_gated2 * (lf.z_esm_cal - lf.z_mkmp_cal) + lf.b_gated2
    lf["p_learned"] = sigmoid(z)
    gt = lf.copy()
    gt["plddt"] = gt["name"].map(plddt_map)
    gt["identity"] = gt["name"].map(hom)
    test_name_to_y = {data_list[i].name: int(data_list[i].y) for i in idx["test"]}
    gt["y"] = gt["name"].map(test_name_to_y)
    gt = gt.dropna(subset=["y", "identity", "plddt"])
    gt["p_explicit"] = apply_explicit(theta, gt.z_mkmp_cal, gt.z_esm_cal, gt.plddt, gt.identity,
                                      val.plddt.mean(), val.identity.mean(), val.plddt.std(), val.identity.std())
    gt["p_mkmp"] = sigmoid(gt.z_mkmp_cal)
    gt["p_esm"] = sigmoid(gt.z_esm_cal)

    def report(pcol):
        best = np.maximum(mcc(gt.y, gt.p_mkmp), mcc(gt.y, gt.p_esm))
        d = {}
        d["overall_mcc"] = mcc(gt.y, gt[pcol])
        d["overall_auc"] = roc_auc_score(gt.y, gt[pcol])
        d["dMCC"] = d["overall_mcc"] - best
        for plo, phi, plab in [(-1, 70, "<70"), (70, 90, "70-90"), (90, 101, ">90")]:
            for ilo, ihi, ilab in [(-0.001, 0.20, "<0.20"), (0.20, 0.40, "0.20-0.40"), (0.40, 1.01, ">0.40")]:
                cell = gt[(gt.plddt >= plo) & (gt.plddt < phi) & (gt.identity >= ilo) & (gt.identity < ihi)]
                if len(cell) < 20:
                    continue
                bm = max(mcc(cell.y, cell.p_mkmp), mcc(cell.y, cell.p_esm))
                d[f"{plab}_{ilab}_dMCC"] = mcc(cell.y, cell[pcol]) - bm
        return d

    learned = report("p_learned")
    explicit = report("p_explicit")
    print("\n== learned gate ==")
    for k, v in learned.items():
        print(f"  {k}: {round(v, 4)}")
    print("\n== explicit gate ==")
    for k, v in explicit.items():
        print(f"  {k}: {round(v, 4)}")

    import pandas as pd
    summ = pd.DataFrame({"metric": list(learned.keys()),
                         "learned_gate": [learned[k] for k in learned.keys()],
                         "explicit_gate": [explicit[k] for k in learned.keys()]})
    out = os.path.join(RESDIR, "e24_explicit_gate_valid_summary.csv")
    summ.to_csv(out, index=False)
    print("\nsaved:", out)


if __name__ == "__main__":
    main()
