"""run_nmr.py

Entry script: sarcopenia prediction from NMR metabolites combined with
metagenomics CLR taxa and clinical covariates. The NMR loader threads three
scaled feature groups (continuous covariates, CLR taxa, NMR metabolites)
through the shared engine.

    python run_nmr.py
"""
import pandas as pd

import config
from loaders import load_nmr
from ml_engine import run_model_set
from feature_importance import (
    lasso_coef_summary, boruta_rf_summary, fcnn_ig_summary,
)

cfg = config.NMR
out = cfg["output_dir"]
out.mkdir(parents=True, exist_ok=True)

meta_filtered = pd.read_csv(cfg["input_dir"] / cfg["meta_filtered_csv"])
meta_df = pd.read_csv(cfg["input_dir"] / cfg["meta_df_csv"])
nmr_df = pd.read_csv(cfg["nmr_dir"] / cfg["nmr_csv"])

# Attach File_ID to the NMR table (loader indexes metabolites by File_ID).
nmr_df = pd.merge(nmr_df, meta_df[["record_id", "File_ID"]],
                  on="record_id", how="left")

X, y, feature_groups, binary_cols = load_nmr(
    meta_filtered, meta_df, nmr_df,
    start_col=cfg["nmr_start_col"], end_col=cfg["nmr_end_col"],
    drop_metabolites=cfg["nmr_drop_metabolites"],
    log_transform=cfg["log_transform_nmr"])

# For the *_only run, the assay groups are taxa + NMR (covariates excluded);
# run_model_set infers the first non-cont group name for file labeling.
run_model_set(X, y, feature_groups, binary_cols, out,
              feature_group_name="clr_taxa",
              outer_seed=config.OUTER_SEED, inner_seed=config.INNER_SEED,
              scoring_metric=config.SCORING_METRIC)

# --- Feature importance from the saved CV results (full models) ---
lasso_coef_summary(
    pd.read_pickle(out / "cv_results_lasso_logreg_full.pkl"),
    outfile=out / "lasso_L1_stratifiedkfold_nmr_i111_o222_coef_summary.csv")
boruta_rf_summary(
    pd.read_pickle(out / "cv_results_random_forest_full.pkl"),
    X, y, outfile=out / "rf_feature_importances_Boruta_nmr_5stratkfold.csv")
fcnn_ig_summary(
    pd.read_pickle(out / "cv_results_fcnn_full.pkl"),
    X.astype("float32"), outfile=out / "fcnn_feature_importances_IG_nmr.csv")
print("NMR ML run complete.")
