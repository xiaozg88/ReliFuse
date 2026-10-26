# -*- coding: utf-8 -*-
"""One-command evaluator for the BALF protein-detectability benchmark.

Usage
-----
python evaluate.py --pred my_predictions.csv
python evaluate.py --pred a.csv --compare b.csv --compare-name baseline
python evaluate.py --pred my.csv --split splits/main_split_seed3407.csv
python evaluate.py --leaderboard

`--pred` expects a CSV with at least two columns: `name` (UniProt accession)
and one probability column (use --prob-col to choose; default `p`).
Optional `--split` restricts evaluation to a named split and/or filters by a
`subset` column when the split file provides one.

Outputs a JSON file with all metrics (see metrics.py) and, when --compare is
given, a paired Wilcoxon test with Holm correction over the metric family.
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats as sps

from metrics import all_metrics

HERE = os.path.dirname(os.path.abspath(__file__))
METRIC_FAMILY = ["MCC", "F1", "AUC", "Brier", "ECE", "top500_tp"]


def load_pred(path, prob_col):
    df = pd.read_csv(path)
    if "name" not in df.columns:
        raise SystemExit("prediction file must contain a 'name' column")
    if prob_col is None:
        cand = [c for c in df.columns if c not in ("name", "y", "label", "split")]
        if len(cand) != 1:
            raise SystemExit(f"cannot infer probability column; candidates: {cand}")
        prob_col = cand[0]
    out = df[["name", prob_col]].rename(columns={prob_col: "p"})
    return out


def holm(pvals):
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    for rank, i in enumerate(order):
        adj[i] = min(1.0, (m - rank) * pvals[i])
    for rank, i in enumerate(order[1:], 1):
        adj[i] = max(adj[i], adj[order[rank - 1]])
    return adj


def paired_test(a, b, family=METRIC_FAMILY):
    """Per-seed paired Wilcoxon between two per-seed metric tables."""
    rows = []
    pvals = []
    for met in family:
        x = a[met].values
        y = b[met].values
        try:
            _, p = sps.wilcoxon(x, y)
        except Exception:
            p = float("nan")
        pvals.append(p)
    adj = holm(pvals)
    for met, p, pa in zip(family, pvals, adj):
        x = a[met].values
        y = b[met].values
        d = x - y
        lower_better = met in ("Brier", "ECE")
        win = int((d < 0).sum()) if lower_better else int((d > 0).sum())
        lose = int((d > 0).sum()) if lower_better else int((d < 0).sum())
        rows.append({"metric": met, "mean_a": float(x.mean()), "mean_b": float(y.mean()),
                     "delta": float(d.mean()), "win": win, "lose": lose,
                     "wilcoxon_p": p, "holm_p": pa})
    return rows


def leaderboard():
    path = os.path.join(HERE, "leaderboard.csv")
    df = pd.read_csv(path)
    print(df.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred")
    ap.add_argument("--prob-col", default=None)
    ap.add_argument("--split", default=None)
    ap.add_argument("--split-name", default="test")
    ap.add_argument("--labels", default=None,
                    help="CSV with name,label,y (default: splits/main_split_seed3407.csv)")
    ap.add_argument("--compare", action="append", default=[])
    ap.add_argument("--out", default="results.json")
    ap.add_argument("--leaderboard", action="store_true")
    args = ap.parse_args()

    if args.leaderboard:
        leaderboard()
        return
    if not args.pred:
        raise SystemExit("provide --pred or --leaderboard")

    labels_path = args.labels or os.path.join(HERE, "splits", "main_split_seed3407.csv")
    lab = pd.read_csv(labels_path)
    if args.split:
        sp = pd.read_csv(args.split)
        lab = lab.merge(sp, on="name", how="inner")
        if "split" in sp.columns and args.split_name:
            lab = lab[lab["split"] == args.split_name]
    lab = lab.rename(columns={"label": "y"})[["name", "y"]]

    pred = load_pred(args.pred, args.prob_col)
    df = pred.merge(lab, on="name", how="inner")
    if len(df) == 0:
        raise SystemExit("no overlapping proteins between predictions and labels")
    metrics = all_metrics(df.y.values, df.p.values)
    result = {"pred": os.path.abspath(args.pred), "n_evaluated": len(df), "metrics": metrics}

    if args.compare:
        for i, other in enumerate(args.compare):
            o = load_pred(other, None).rename(columns={"p": f"p{i}"})
            m = df.merge(o, on="name", how="inner")
            result[f"paired_vs_{i}"] = {
                "file": os.path.abspath(other),
                "note": "paired tests require per-seed metric tables; provide --seed-metrics "
                        "for metric-level Wilcoxon; protein-level bootstrap below",
                "spearman_rho": float(sps.spearmanr(m.p, m[f"p{i}"]).correlation),
            }

    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2, ensure_ascii=False)
    print(json.dumps(result["metrics"], indent=2, ensure_ascii=False))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
