"""
extract_fold_auc.py

Recovers correct per-outer-fold ROC AUC from ALREADY-TRAINED cv_results
pickles, without retraining anything -- just inference on each fold's held-out
samples using the estimator already stored in cv_results["estimator"].

This sidesteps the ml_engine.py bug where cv_results["test_score"] was saved
using the default scorer (accuracy) instead of roc_auc: the ROC-curve PDFs
were always computing the correct AUC via RocCurveDisplay/roc_auc_score, this
script just does the same computation and saves it to a CSV instead of only
drawing it on a plot.

Requires the SAME X, y used at training time -- reconstructed here via
loaders.load_metagenomics(), which only re-reads the metadata/abundance CSVs
and re-does the (cheap) CLR/log transform; it does not refit any model.

Usage:
  python extract_fold_auc.py <results_dir> \
      --meta-filtered-csv meta_filtered_06242026.csv \
      --meta-df-csv meta_df_FullSaMu_06242026.csv \
      --input-dir /data/local/jy1008/SaMu/results/latest/metagenomics_R \
      --models lasso_logreg random_forest fcnn \
      --feature-group clr_taxa \
      --out fold_auc.csv
"""
import argparse
from pathlib import Path

import pandas as pd
from sklearn.metrics import roc_auc_score

from loaders import load_metagenomics

MODEL_DISPLAY = {
    "lasso_logreg": "Lasso logistic\nregression",
    "random_forest": "Random forest",
    "fcnn": "FCNN",
    "mlp": "MLP",
    "adaboost": "AdaBoost",
}

# Models trained on CLR-transformed features vs. log-transformed (matches
# run_metagenomics.py: lasso_logreg/fcnn use X_clr, random_forest uses X_log).
CLR_MODELS = {"lasso_logreg", "fcnn", "mlp", "adaboost"}


def fold_aucs_for_pickle(path, X, y):
    """Recompute per-fold AUC from a cv_results pickle's stored estimators,
    the same way RocCurveDisplay.from_cv_results does internally."""
    cv_results = pd.read_pickle(path)
    aucs = []
    for fold_i, (estimator, test_idx) in enumerate(
            zip(cv_results["estimator"], cv_results["indices"]["test"]), 1):
        X_test = X.iloc[test_idx]
        y_test = y.iloc[test_idx]
        y_score = estimator.predict_proba(X_test)[:, 1]
        aucs.append((fold_i, roc_auc_score(y_test, y_score)))
    return aucs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("--input-dir", required=True,
                     help="Directory containing meta_filtered/meta_df CSVs "
                          "(same as config.METAGENOMICS['input_dir'])")
    ap.add_argument("--meta-filtered-csv", required=True)
    ap.add_argument("--meta-df-csv", required=True)
    ap.add_argument("--models", nargs="+",
                     default=["lasso_logreg", "random_forest", "fcnn"])
    ap.add_argument("--feature-group", default="clr_taxa")
    ap.add_argument("--out", default="fold_auc.csv")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    input_dir = Path(args.input_dir)

    meta_filtered = pd.read_csv(input_dir / args.meta_filtered_csv)
    meta_df = pd.read_csv(input_dir / args.meta_df_csv)

    # Reconstruct both feature matrices exactly as run_metagenomics.py did --
    # inference-only, no fitting.
    X_clr, y, *_ = load_metagenomics(meta_filtered, meta_df, clr_transform=True)
    X_log, _, *_ = load_metagenomics(meta_filtered, meta_df, clr_transform=False)

    subset_display = {
        "full": "All features",
        "clinical_only": "Covariates only",
        f"{args.feature_group}_only": "Microbiome only",
    }

    rows = []
    for model in args.models:
        if model == "fcnn":
            X = X_clr.astype("float32")
        elif model in CLR_MODELS:
            X = X_clr
        else:
            X = X_log
        for subset_key, subset_label in subset_display.items():
            path = results_dir / f"cv_results_{model}_{subset_key}.pkl"
            if not path.exists():
                print(f"Skipping missing file: {path.name}")
                continue
            for fold_i, auc in fold_aucs_for_pickle(path, X, y):
                rows.append({
                    "model": MODEL_DISPLAY.get(model, model),
                    "subset": subset_label,
                    "fold": fold_i,
                    "score": auc,
                })

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(f"Saved {len(df)} fold-level AUC rows to {args.out}")
    print(df.groupby(["model", "subset"])["score"].agg(["mean", "std"]))


if __name__ == "__main__":
    main()