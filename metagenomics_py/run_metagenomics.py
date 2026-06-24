"""run_metagenomics.py

Entry script: sarcopenia prediction from MetaPhlAn species relative abundances
plus clinical covariates. Thin wrapper over the shared engine -- it loads the
metagenomics data, runs the model sweep, and extracts feature importances.

    python run_metagenomics.py
"""
import pandas as pd

import config
from loaders import load_metagenomics
from ml_engine import run_model_set
from feature_importance import (
    lasso_coef_summary, boruta_rf_summary, fcnn_ig_summary,
    merge_importance_tables,
)

cfg = config.METAGENOMICS
out = cfg["output_dir"]
out.mkdir(parents=True, exist_ok=True)

meta_filtered = pd.read_csv(cfg["input_dir"] / cfg["meta_filtered_csv"])
meta_df = pd.read_csv(cfg["input_dir"] / cfg["meta_df_csv"])

# --- Lasso + FC_NN use CLR-transformed taxa; random forest uses log-abundance.
X_clr, y, groups_clr, binary_cols = load_metagenomics(
    meta_filtered, meta_df, clr_transform=True)
X_log, _, groups_log, _ = load_metagenomics(
    meta_filtered, meta_df, clr_transform=False)

run_model_set(X_clr, y, groups_clr, binary_cols, out,
              models=("lasso_logreg", "fcnn"), feature_group_name="clr_taxa",
              outer_seed=config.OUTER_SEED, inner_seed=config.INNER_SEED,
              scoring_metric=config.SCORING_METRIC)
run_model_set(X_log, y, groups_log, binary_cols, out,
              models=("random_forest",), feature_group_name="clr_taxa",
              outer_seed=config.OUTER_SEED, inner_seed=config.INNER_SEED,
              scoring_metric=config.SCORING_METRIC)

# --- Feature importance from the saved CV results (full models) ---
linreg_res = lasso_coef_summary(
    pd.read_pickle(out / "cv_results_lasso_logreg_full.pkl"),
    outfile=out / "lasso_L1_stratifiedkfold_i111_o222_coef_summary.csv")

rf_res = boruta_rf_summary(
    pd.read_pickle(out / "cv_results_random_forest_full.pkl"),
    X_log, y, outfile=out / "rf_feature_importances_Boruta_5stratkfold.csv")

fcnn_res = fcnn_ig_summary(
    pd.read_pickle(out / "cv_results_fcnn_full.pkl"),
    X_clr.astype("float32"),
    outfile=out / "fcnn_feature_importances_IG.csv")

# --- Merge with DESeq2 into one feature table ---
deseq2_res = pd.read_csv(cfg["input_dir"] / cfg["deseq2_csv"])
merge_importance_tables(deseq2_res, linreg_res, rf_res, fcnn_res,
                        outfile=out / "merged_feature_importance.csv")
print("Metagenomics ML run complete.")
