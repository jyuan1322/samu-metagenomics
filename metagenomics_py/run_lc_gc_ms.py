"""run_lc_gc_ms.py

Entry script: sarcopenia prediction from LC-MS / GC-MS features (DEP2
VSN-imputed) plus clinical covariates. Set the assay via
config.LC_GC_MS["experiment_name"] ("GC_MS", "LC_MS_neg", "LC_MS_pos").

    python run_lc_gc_ms.py
"""
import pandas as pd

import config
from loaders import load_lc_gc_ms
from ml_engine import run_model_set
from feature_importance import (
    lasso_coef_summary, boruta_rf_summary, fcnn_ig_summary,
)

cfg = config.LC_GC_MS
exp = cfg["experiment_name"]
out = cfg["output_dir"]
out.mkdir(parents=True, exist_ok=True)

matrix_csv = cfg["matrix_csv_tmpl"].format(experiment_name=exp)
meta_csv = cfg["meta_csv_tmpl"].format(experiment_name=exp)
feature_df = pd.read_csv(cfg["input_dir"] / matrix_csv)
meta_df = pd.read_csv(cfg["input_dir"] / meta_csv)

X, y, feature_groups, binary_cols = load_lc_gc_ms(feature_df, meta_df)

run_model_set(X, y, feature_groups, binary_cols, out,
              feature_group_name="gcms",
              outer_seed=config.OUTER_SEED, inner_seed=config.INNER_SEED,
              scoring_metric=config.SCORING_METRIC)

# --- Feature importance from the saved CV results (full models) ---
lasso_coef_summary(
    pd.read_pickle(out / "cv_results_lasso_logreg_full.pkl"),
    outfile=out / f"lasso_coef_summary_{exp}.csv")
boruta_rf_summary(
    pd.read_pickle(out / "cv_results_random_forest_full.pkl"),
    X, y, outfile=out / f"rf_boruta_summary_{exp}.csv")
fcnn_ig_summary(
    pd.read_pickle(out / "cv_results_fcnn_full.pkl"),
    X.astype("float32"), outfile=out / f"fcnn_feature_importances_IG_{exp}.csv")
print(f"LC/GC-MS ({exp}) ML run complete.")
