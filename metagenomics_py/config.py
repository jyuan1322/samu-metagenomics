"""config.py

Paths and shared run settings for the SaMu ML entry scripts. Edit the paths
here per study; the entry scripts and engine read from this module so no paths
are hardcoded in the run logic.
"""
from pathlib import Path

# ---------------------------------------------------------------------------
# Nested-CV seeds (outer/inner StratifiedKFold). Kept here so every assay runs
# the same split scheme; change in one place to re-run under different seeds.
# ---------------------------------------------------------------------------
OUTER_SEED = 222
INNER_SEED = 111
SCORING_METRIC = "roc_auc"

# ---------------------------------------------------------------------------
# Metagenomics
# ---------------------------------------------------------------------------
METAGENOMICS = {
    "input_dir": Path("/data/local/jy1008/SaMu/results/latest/metagenomics_R"),
    "output_dir": Path("/data/local/jy1008/SaMu/results/latest/metagenomics_ml"),
    "meta_filtered_csv": "meta_filtered_06242026.csv",
    "meta_df_csv": "meta_df_FullSaMu_06242026.csv",
    "deseq2_csv": "deseq2_results_06242026.csv",
}

# ---------------------------------------------------------------------------
# Mass-spec proteomics
# ---------------------------------------------------------------------------
MASS_SPEC = {
    "input_dir": Path("/data/local/jy1008/SaMu/proteomics/results"),
    "meta_dir": Path("/data/local/jy1008/SaMu/results/02102026"),
    "output_dir": Path("/data/local/jy1008/SaMu/results/latest/proteomics_mass_spec/results"),
    "matrix_csv": "proteomics_dep2_vsn_imputed_matrix_p15filt.csv",
    "meta_df_csv": "meta_df_FullSaMu.csv",
}

# ---------------------------------------------------------------------------
# NMR (combined with metagenomics taxa)
# ---------------------------------------------------------------------------
NMR = {
    "input_dir": Path("/data/local/jy1008/SaMu/results/02102026"),
    "nmr_dir": Path("/data/local/jy1008/SaMu/metadata"),
    "output_dir": Path("/data/local/jy1008/SaMu/results/latest/nmr_ml"),
    "meta_filtered_csv": "meta_filtered.csv",
    "meta_df_csv": "meta_df_FullSaMu.csv",
    "nmr_csv": "20251027_SaMu_NMR_ureum_creat_v10.csv",
    "nmr_start_col": "5-Aminopentanoate",
    "nmr_end_col": "Valine",
    "nmr_drop_metabolites": ("Glucose", "Galactose"),
    "log_transform_nmr": True,
}

# ---------------------------------------------------------------------------
# LC-MS / GC-MS. Set experiment_name to the assay; paths derive from it.
#   experiment_name in {"GC_MS", "LC_MS_neg", "LC_MS_pos"}
# ---------------------------------------------------------------------------
LC_GC_MS = {
    "experiment_name": "GC_MS",
    "input_dir": Path("/data/local/jy1008/SaMu/results/latest/proteomics_GC-MS"),
    "output_dir": Path("/data/local/jy1008/SaMu/results/latest/proteomics_GC-MS/ml_models"),
    # matrix/metadata filenames are f-string templated on experiment_name:
    "matrix_csv_tmpl": "proteomics_{experiment_name}_vsn_imputed_matrix.csv",
    "meta_csv_tmpl": "proteomics_{experiment_name}_vsn_imputed_metadata.csv",
}
