"""run_mass_spec.py

Entry script: sarcopenia prediction from mass-spec proteomics (DEP2 VSN-imputed
intensities) plus clinical covariates.

    python run_mass_spec.py
"""
import pandas as pd

import config
from loaders import load_mass_spec
from ml_engine import run_model_set
from feature_importance import (
    lasso_coef_summary, boruta_rf_summary, fcnn_ig_summary,
)

cfg = config.MASS_SPEC
out = cfg["output_dir"]
out.mkdir(parents=True, exist_ok=True)

data_filtered = pd.read_csv(cfg["input_dir"] / cfg["matrix_csv"])
meta_df = pd.read_csv(cfg["meta_dir"] / cfg["meta_df_csv"])

X, y, feature_groups, binary_cols = load_mass_spec(data_filtered, meta_df)

run_model_set(X, y, feature_groups, binary_cols, out,
              feature_group_name="proteins",
              outer_seed=config.OUTER_SEED, inner_seed=config.INNER_SEED,
              scoring_metric=config.SCORING_METRIC)

# --- Feature importance from the saved CV results (full models) ---
lasso_coef_summary(
    pd.read_pickle(out / "cv_results_lasso_logreg_full.pkl"),
    outfile=out / "lasso_L1_stratifiedkfold_i111_o222_coef_summary.csv")
boruta_rf_summary(
    pd.read_pickle(out / "cv_results_random_forest_full.pkl"),
    X, y, outfile=out / "rf_feature_importances_Boruta_5stratkfold.csv")
fcnn_ig_summary(
    pd.read_pickle(out / "cv_results_fcnn_full.pkl"),
    X.astype("float32"), outfile=out / "fcnn_feature_importances_IG.csv")
print("Mass-spec ML run complete.")
