# -*- coding: utf-8 -*-
"""External-BALF recalibration: scale drift vs feature failure (Hmmier et al. 2017, BBA Clin 7:97-104, PMC5357681).

Question posed in the v51 work log: is the cross-cohort negative transfer of
the fusion model caused by *structure-feature failure* or by *probability
scale drift*?  This experiment separates the two on the protein-disjoint
Hmmier 2017 subset (n = 270; 91 positives / 179 negatives; the same cohort as
Supplementary Table 6 / main-text Table XII).

Protocol
--------
500 stratified random 50/50 calibration/evaluation splits (seed 3407).  For
each split, recalibration operators are fitted on the calibration half only
and evaluated once on the evaluation half:

- zero-shot: raw sequence / structure / ReliFuse probabilities (no fitting);
- temperature scaling (1 parameter) on each pathway logit;
- Platt scaling (2 parameters) on each pathway logit;
- isotonic regression on the sequence probability;
- recalibrated fusion: per-pathway Platt probabilities averaged;
- two-pathway logistic stacking on the raw calibrated logits.

Interpretation
--------------
- Monotone recalibration cannot change a single score's ranking, so
  structure-only AUC after recalibration isolates *feature failure*.
- Per-pathway recalibration before fusion changes the relative influence of
  the two pathways, so the fused AUC after recalibration isolates how much of
  the negative transfer is a *probability-scale / saturation* artifact.

Output: results/e26_external_recalibration.csv
"""

import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import wilcoxon
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, matthews_corrcoef, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
RESDIR = os.path.join(ROOT, "results")
SRC = os.path.join(ROOT, "data", "bal_ext_predictions.csv")
SPLIT = os.path.join(ROOT, "data", "splits", "main_split_seed3407.csv")

N_REPS = 500
SEED = 3407
EPS = 1e-6


def load_disjoint():
    ext = pd.read_csv(SRC)
    bench = set(pd.read_csv(SPLIT)["name"].astype(str))
    d = ext[~ext["accession"].astype(str).isin(bench)].reset_index(drop=True)
    assert d["detect"].value_counts().to_dict() == {0: 179, 1: 91}
    return d


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -60.0, 60.0)))


def ece(y, p, nbins=10):
    e = 0.0
    for k in range(nbins):
        lo, hi = k / nbins, (k + 1) / nbins
        m = (p >= lo) & ((p < hi) if k < nbins - 1 else (p <= hi))
        if m.sum():
            e += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(e)


def metrics(y, p):
    p = np.clip(p, EPS, 1 - EPS)
    return {"AUC": roc_auc_score(y, p),
            "MCC": matthews_corrcoef(y, (p >= 0.5).astype(int)),
            "Brier": brier_score_loss(y, p),
            "ECE": ece(y, p)}


def fit_temperature(z, y):
    """1-parameter temperature scaling on a logit (NLL, log-parameterized)."""
    def nll(log_t):
        p = sigmoid(z / np.exp(log_t))
        return -np.mean(y * np.log(p + EPS) + (1 - y) * np.log(1 - p + EPS))
    res = minimize_scalar(nll, bounds=(-6.0, 12.0), method="bounded")
    return np.exp(res.x)


def fit_platt(z, y):
    m = LogisticRegression(C=1e6, max_iter=2000)
    m.fit(z.reshape(-1, 1), y)
    return lambda zz: m.predict_proba(zz.reshape(-1, 1))[:, 1]


def main():
    d = load_disjoint()
    y = d["detect"].values.astype(int)
    zs = d["z_mkmp_cal"].values.astype(float)
    zq = d["z_esm_cal"].values.astype(float)
    ps_raw = d["p_mkmp"].values.astype(float)
    pq_raw = d["p_esm"].values.astype(float)
    pf_raw = d["p_lgf"].values.astype(float)

    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    rng = np.random.default_rng(SEED)
    rows = []
    for rep in range(N_REPS):
        cal = np.concatenate([rng.choice(pos, 46, replace=False),
                              rng.choice(neg, 90, replace=False)])
        ev = np.setdiff1d(np.arange(len(y)), cal)

        # --- single-score recalibration operators (fitted on cal) ---
        T_s = fit_temperature(zs[cal], y[cal])
        T_q = fit_temperature(zq[cal], y[cal])
        ms = LogisticRegression(C=1e6, max_iter=2000).fit(zs[cal].reshape(-1, 1), y[cal])
        mq = LogisticRegression(C=1e6, max_iter=2000).fit(zq[cal].reshape(-1, 1), y[cal])
        p_platt_s = ms.predict_proba(zs[ev].reshape(-1, 1))[:, 1]
        p_platt_q = mq.predict_proba(zq[ev].reshape(-1, 1))[:, 1]
        p_temp_s = sigmoid(zs[ev] / T_s)
        p_temp_q = sigmoid(zq[ev] / T_q)
        iso = IsotonicRegression(out_of_bounds="clip").fit(pq_raw[cal], y[cal])
        p_iso_q = iso.predict(pq_raw[ev])

        # --- fusion after pathway-level recalibration ---
        p_recal_mean = 0.5 * (p_platt_s + p_platt_q)
        stack = LogisticRegression(C=1e6, max_iter=2000).fit(
            np.column_stack([zs[cal], zq[cal]]), y[cal])
        p_stack = stack.predict_proba(np.column_stack([zs[ev], zq[ev]]))[:, 1]

        entries = {
            "zero_sequence": pq_raw[ev],
            "zero_structure": ps_raw[ev],
            "zero_fusion": pf_raw[ev],
            "temp_structure": p_temp_s,
            "temp_sequence": p_temp_q,
            "platt_structure": p_platt_s,
            "platt_sequence": p_platt_q,
            "isotonic_sequence": p_iso_q,
            "recalibrated_mean_fusion": p_recal_mean,
            "stacking_2path": p_stack,
        }
        for name, p in entries.items():
            rows.append({"rep": rep, "method": name, **metrics(y[ev], p)})

    r = pd.DataFrame(rows)

    # aggregate
    agg = (r.groupby("method")[["AUC", "MCC", "Brier", "ECE"]]
             .agg(["mean", "std"]).round(4))
    agg.columns = [f"{m}_{s}" for m, s in agg.columns]
    agg = agg.reset_index()

    # paired comparisons vs zero-shot fusion and zero-shot sequence
    piv = {m: g.sort_values("rep").reset_index(drop=True)
           for m, g in r.groupby("method")}
    pair_rows = []
    for metric in ["AUC", "MCC", "Brier", "ECE"]:
        for ref in ["zero_fusion", "zero_sequence"]:
            for m in ["platt_sequence", "isotonic_sequence",
                      "recalibrated_mean_fusion", "stacking_2path",
                      "platt_structure"]:
                delta = piv[m][metric] - piv[ref][metric]
                if float(np.abs(delta).max()) == 0.0:
                    p = 1.0
                else:
                    p = float(wilcoxon(piv[m][metric], piv[ref][metric]).pvalue)
                pair_rows.append({"method": m, "reference": ref, "metric": metric,
                                  "mean_delta": round(float(delta.mean()), 4),
                                  "win/loss/tie": f"{int((delta>0).sum())}/{int((delta<0).sum())}/{int((delta==0).sum())}",
                                  "wilcoxon_p": float(p)})

    os.makedirs(RESDIR, exist_ok=True)
    out_main = os.path.join(RESDIR, "e26_external_recalibration.csv")
    agg.to_csv(out_main, index=False)
    out_pair = os.path.join(RESDIR, "e26_external_recalibration_paired.csv")
    pd.DataFrame(pair_rows).to_csv(out_pair, index=False)
    print(f"saved {out_main}")
    print(f"saved {out_pair}")
    print(agg.to_string(index=False))


if __name__ == "__main__":
    main()

