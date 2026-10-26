# -*- coding: utf-8 -*-
"""Unified evaluation metrics for the BALF protein-detectability benchmark.

All metrics are computed at the protein level on a binary detection label.
Calibration (Brier, ECE) and decision-utility (top-K, fixed-specificity
sensitivity, DCA net benefit) metrics are first-class citizens: a model that
ranks well but is miscalibrated is penalised on the utility side.
"""

import numpy as np
from sklearn.metrics import (accuracy_score, average_precision_score,
                             brier_score_loss, f1_score, matthews_corrcoef,
                             precision_score, recall_score, roc_auc_score)


def ece(y, p, n_bins=10):
    """Expected calibration error, 10 equal-width bins, last bin closed."""
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    out = 0.0
    for i in range(n_bins):
        if i < n_bins - 1:
            m = (p >= edges[i]) & (p < edges[i + 1])
        else:
            m = (p >= edges[i]) & (p <= edges[i + 1])
        if m.sum():
            out += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(out)


def topk_metrics(y, p, ks=(50, 100, 200, 500, 1000)):
    y = np.asarray(y)
    p = np.asarray(p)
    order = np.argsort(-p)
    out = {}
    for k in ks:
        kk = min(k, len(y))
        tp = int(y[order[:kk]].sum())
        out[f"top{k}_tp"] = tp
        out[f"top{k}_precision"] = tp / kk
        out[f"top{k}_recall"] = tp / max(int(y.sum()), 1)
    return out


def sensitivity_at_specificity(y, p, target=0.99):
    y = np.asarray(y)
    p = np.asarray(p)
    neg = p[y == 0]
    thr = np.quantile(neg, target) if len(neg) else 0.5
    pred = (p >= thr).astype(int)
    return float(recall_score(y, pred, zero_division=0))


def precision_at_recall(y, p, target=0.9):
    y = np.asarray(y)
    p = np.asarray(p)
    order = np.argsort(-p)
    ys = y[order]
    cum = np.cumsum(ys)
    total = max(int(y.sum()), 1)
    idx = np.searchsorted(cum, target * total)
    idx = min(int(idx), len(ys) - 1)
    return float(ys[:idx + 1].mean())


def dca_net_benefit(y, p, thresholds=(0.2, 0.3, 0.4, 0.5)):
    """Net benefit = TP/N - FP/N * pt/(1-pt)."""
    y = np.asarray(y)
    p = np.asarray(p)
    n = len(y)
    out = {}
    for t in thresholds:
        pred = (p >= t).astype(int)
        tp = int(((pred == 1) & (y == 1)).sum())
        fp = int(((pred == 1) & (y == 0)).sum())
        out[f"nb@{t}"] = tp / n - (fp / n) * (t / (1 - t))
    return out


def all_metrics(y, p, ks=(50, 100, 200, 500, 1000)):
    y = np.asarray(y)
    p = np.asarray(p)
    pred = (p >= 0.5).astype(int)
    m = {
        "n": int(len(y)),
        "n_pos": int(y.sum()),
        "AUC": float(roc_auc_score(y, p)),
        "PR_AUC": float(average_precision_score(y, p)),
        "Accuracy": float(accuracy_score(y, pred)),
        "MCC": float(matthews_corrcoef(y, pred)),
        "F1": float(f1_score(y, pred)),
        "Precision": float(precision_score(y, pred, zero_division=0)),
        "Recall": float(recall_score(y, pred, zero_division=0)),
        "Brier": float(brier_score_loss(y, p)),
        "ECE": ece(y, p),
        "Sens@Spec0.99": sensitivity_at_specificity(y, p, 0.99),
        "Prec@Recall0.9": precision_at_recall(y, p, 0.9),
    }
    m.update(topk_metrics(y, p, ks))
    m.update(dca_net_benefit(y, p))
    return m
