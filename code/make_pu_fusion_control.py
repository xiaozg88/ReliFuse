# -*- coding: utf-8 -*-
"""Positive-unlabeled (PU) fusion-head control on cached pathway logits.

The v51 work log requested a method-level PU comparison instead of purely
evaluation-level PU sensitivity.  Using the 20-seed validation/test calibrated
logit cache (structure z_mkmp_cal, sequence z_esm_cal), this script fits
alternative PU-risk fusion heads on the VALIDATION split only and evaluates
once on the TEST split, per seed:

  - bce       : observed-label BCE (the current ReliFuse training convention,
                instantiated here as a linear logistic head for comparability);
  - upu       : unbiased PU risk estimator (du Plessis et al., 2015),
                R = pi*E_P[softplus(-f)] + E_U[softplus(f)] - pi*E_P[softplus(f)];
  - nnpu      : non-negative PU estimator (Kiryo et al., 2017), same risk with
                the correction term clipped at zero;
  - pu_bagging: bagged BCE heads (Mordelet & Vert, 2014), each trained on all
                labeled positives plus a fresh random half of the unlabeled set,
                10 members, averaged probabilities;
  - conf2step : two-step positive-confidence weighting - a stage-1 BCE head
                provides soft targets for the unlabeled data, and a stage-2 head
                is trained with those soft targets.

Class-prior sensitivity for upu/nnpu: pi in {0.298 (observed validation
positive rate), 0.40, 0.50}.

Honest evaluation note: the TEST labels are themselves PU (non-detection is
not a biological negative), so this comparison quantifies how PU-corrected
training objectives behave under the paper's standard observational-label
evaluation; it is not a clean-label benchmark.

Output: results/e27_pu_fusion_control.csv (per-seed), 
        results/e27_pu_fusion_control_summary.csv (mean +/- SD).
"""

import os

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.metrics import brier_score_loss, matthews_corrcoef, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
RESDIR = os.path.join(ROOT, "results")
CACHE = os.path.join(RESDIR, "val_logits_cache.csv")
SPLIT = os.path.join(ROOT, "data", "splits", "main_split_seed3407.csv")

EPOCHS = 1000
LR = 0.05
EPS = 1e-6


def ece(y, p, nbins=10):
    e = 0.0
    for k in range(nbins):
        lo, hi = k / nbins, (k + 1) / nbins
        m = (p >= lo) & ((p < hi) if k < nbins - 1 else (p <= hi))
        if m.sum():
            e += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(e)


def sigmoid(z):
    z = np.clip(z, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-z))


def fit_head(Xv, yv, loss, pi=0.298, seed=3407, soft_targets=None, epochs=EPOCHS):
    """Full-batch Adam on a linear logistic head; returns (w, b)."""
    torch.manual_seed(seed)
    w = torch.zeros(Xv.shape[1], requires_grad=True)
    b = torch.zeros(1, requires_grad=True)
    X = torch.tensor(Xv, dtype=torch.float32)
    y = torch.tensor(yv, dtype=torch.float32)
    P, U = y == 1, y == 0
    opt = torch.optim.Adam([w, b], lr=LR)
    for _ in range(epochs):
        opt.zero_grad()
        f = X @ w + b
        if loss == "bce":
            t = y if soft_targets is None else torch.tensor(soft_targets, dtype=torch.float32)
            lossv = F.binary_cross_entropy_with_logits(f, t)
        elif loss in ("upu", "nnpu"):
            T1 = F.softplus(-f[P]).mean()   # E_P softplus(-f): P as positive
            T2 = F.softplus(f[U]).mean()    # E_U softplus(+f): U as negative
            T3 = F.softplus(f[P]).mean()    # E_P softplus(+f): correction term
            corr = T2 - pi * T3
            lossv = pi * T1 + (corr if loss == "upu" else torch.clamp(corr, min=0.0))
        lossv.backward()
        torch.nn.utils.clip_grad_norm_([w, b], 10.0)
        opt.step()
    return w.detach().numpy(), float(b.detach().item())


def eval_head(w, b, Xt, yt):
    p = sigmoid(Xt @ w + b)
    return {"AUC": roc_auc_score(yt, p),
            "MCC": matthews_corrcoef(yt, (p >= 0.5).astype(int)),
            "Brier": brier_score_loss(yt, p),
            "ECE": ece(yt, p)}


def main():
    cache = pd.read_csv(CACHE)
    labels = pd.read_csv(SPLIT)[["name", "label"]]
    d = cache.merge(labels, on="name", how="left")
    assert d["label"].notna().all()
    seeds = sorted(d["seed"].unique())

    rows = []
    rng_master = np.random.default_rng(3407)
    for seed in seeds:
        g = d[d.seed == seed]
        val, test = g[g.split == "val"], g[g.split == "test"]
        Xv = val[["z_mkmp_cal", "z_esm_cal"]].values.astype(float)
        yv = val["label"].values.astype(int)
        Xt = test[["z_mkmp_cal", "z_esm_cal"]].values.astype(float)
        yt = test["label"].values.astype(int)
        pi_obs = float(yv.mean())

        configs = [("bce", pi_obs)]
        for pi in (pi_obs, 0.40, 0.50):
            configs.append(("upu", pi))
            configs.append(("nnpu", pi))

        for loss, pi in configs:
            w, b = fit_head(Xv, yv, loss, pi=pi, seed=seed)
            m = eval_head(w, b, Xt, yt)
            rows.append({"seed": seed, "method": loss, "pi": round(pi, 3),
                         "intercept_b": round(b, 3), "w1": round(float(w[0]), 4),
                         "w2": round(float(w[1]), 4), **m})

        # PU bagging: 10 members, all P + fresh random 50% of U as N.
        rng = np.random.default_rng(1000 + int(seed))
        pidx, uidx = np.flatnonzero(yv == 1), np.flatnonzero(yv == 0)
        probs = []
        for m_ in range(10):
            keep = rng.choice(uidx, size=len(uidx) // 2, replace=False)
            idx = np.concatenate([pidx, keep])
            w, b = fit_head(Xv[idx], yv[idx], "bce", seed=seed * 100 + m_)
            probs.append(sigmoid(Xt @ w + b))
        p_bag = np.clip(np.mean(probs, axis=0), EPS, 1 - EPS)
        rows.append({"seed": seed, "method": "pu_bagging", "pi": round(pi_obs, 3),
                     "intercept_b": np.nan, "w1": np.nan, "w2": np.nan,
                     **{k: v for k, v in {
                         "AUC": roc_auc_score(yt, p_bag),
                         "MCC": matthews_corrcoef(yt, (p_bag >= 0.5).astype(int)),
                         "Brier": brier_score_loss(yt, p_bag),
                         "ECE": ece(yt, p_bag)}.items()}})

        # Two-step positive-confidence: stage-1 BCE -> soft targets for U.
        w1v, b1v = fit_head(Xv, yv, "bce", seed=seed)
        soft = yv.astype(float).copy()
        soft[uidx] = sigmoid(Xv[uidx] @ w1v + b1v)
        w2v, b2v = fit_head(Xv, yv, "bce", seed=seed, soft_targets=soft)
        rows.append({"seed": seed, "method": "conf2step", "pi": round(pi_obs, 3),
                     "intercept_b": round(b2v, 3), "w1": round(float(w2v[0]), 4),
                     "w2": round(float(w2v[1]), 4), **eval_head(w2v, b2v, Xt, yt)})

    r = pd.DataFrame(rows)
    os.makedirs(RESDIR, exist_ok=True)
    out = os.path.join(RESDIR, "e27_pu_fusion_control.csv")
    r.to_csv(out, index=False)
    summ = (r.groupby(["method", "pi"])[["AUC", "MCC", "Brier", "ECE", "intercept_b"]]
              .agg(["mean", "std"]).round(4))
    summ.columns = [f"{m}_{s}" for m, s in summ.columns]
    summ = summ.reset_index()
    out2 = os.path.join(RESDIR, "e27_pu_fusion_control_summary.csv")
    summ.to_csv(out2, index=False)
    print(f"saved {out}")
    print(f"saved {out2}")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
