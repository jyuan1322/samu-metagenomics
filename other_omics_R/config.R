# =============================================================================
# config.R — configuration for the SaMu other-omics analyses
# (NMR/quorum metabolomics, DEP2 proteomics/MS, LC-MS feature exploration).
#
# Source at the top of each script, together with the shared metadata utils:
#   source("config.R")
#   source(METADATA_UTILS)
# =============================================================================

# ---------------------------------------------------------------------------
# Metadata utilities
# ---------------------------------------------------------------------------
# Path to metadata_utils.R (local to this repo). The other-omics scripts recode
# from the raw metadata CSV because their assays can include samples not in the
# metagenomics cohort, so they cannot simply reuse the metagenomics output.
METADATA_UTILS <- "metadata_utils.R"

# ---------------------------------------------------------------------------
# Metadata (common to all scripts)
# ---------------------------------------------------------------------------
METADATA_CSV <- "/data/local/jy1008/SaMu/metadata/SaMu_sarcopeniestatus_majorcovariates_v13_16012026.csv"

# Cohort flag column to keep == 1 (NA to skip).
FULLSAMU_COL <- "Full.SaMu"

# Grouping factor used across all analyses (reference level first).
GROUP_VAR    <- "sarc_status_bin"
GROUP_LEVELS <- c("NoSarc", "Sarc")

# Covariate columns retained from the metadata.
META_KEEP_COLS <- c(
  "record_id", "Block", "Full.SaMu", "sarc_status",
  "EWGSOP_strength", "EWGSOP_mass", "EWGSOP_performance",
  "age_def", "sex", "smke", "alco", "nutr_score", "bmi", "stvol"
)

# Set FALSE if the metadata file already stores smke/alco binarized.
RECODE_SMOKE_ALCOHOL <- TRUE

# =============================================================================
# 01_nmr_quorum.R
# =============================================================================
NMR <- list(
  input_dir   = "/data/local/jy1008/SaMu",
  output_dir  = "/data/local/jy1008/SaMu/results/latest/nmr_quorum",

  # Assay tables (relative to input_dir).
  nmr_csv    = "metadata/20251027_SaMu_NMR_ureum_creat_v10.csv",
  quorum_csv = "metadata/20251030_SaMu_QSP_data_v1.csv",

  # NMR metabolite columns start at this 1-based index in nmr_csv.
  nmr_first_col   = 12,
  # Quorum columns are selected by this suffix in quorum_csv.
  quorum_col_regex = "_Quant\\.Prob$",

  # Value standing in for below-detection-limit in NMR data -> 0.
  nmr_below_detection = "UNK_E_1",

  # Metabolites to drop from the Bray-Curtis ordination (high missingness).
  beta_remove_metabs = c("Glucose", "Galactose"),

  epsilon = 1e-6
)

# =============================================================================
# 02_proteomics_dep2.R
# =============================================================================
# loader: "massspec" reads per-sample long CSVs (Protein.Group/Q-value/RawCount,
#         original run_proteomics.R) and runs the contaminant/Q-value QC.
#         "wide" reads sample x feature CSVs that are merged on feature name
#         (original run_proteomics_LC_GC.R); QC is skipped.
PROTEOMICS <- list(
  loader = "massspec",          # "massspec" | "wide"

  parent_dir = "/data/local/jy1008/SaMu/proteomics",
  output_dir = "/data/local/jy1008/SaMu/proteomics/results",
  experiment_name = "proteomics",

  # massspec loader: strip these from each filename to get the sample label.
  sample_name_strip = c("SaMu_sPROT1_rd1_pr1.0_rs0_", "_v1.csv"),
  # massspec QC: keep non-contaminants with Global.Q.Value below this.
  qvalue_cutoff = 0.01,

  # DEP2 missing-value filter: keep features valid in >= `fraction` of samples
  # within at least one condition; `thr` is the allowed missing count.
  filter_thr      = 50,
  filter_fraction = 0.15,
  # Sweep grid for the diagnostic filtering-threshold plot.
  filter_sweep_fractions  = seq(0.1, 0.9, by = 0.05),
  filter_sweep_thresholds = c(10, 25, 50, 100),

  impute_fun = "MinDet",        # DEP2::impute fun

  de_test = "Sarc_vs_NoSarc",   # test_diff contrast (matches GROUP_LEVELS)
  heatmap_top_n = 50            # top proteins by DE p-value in the heatmap
)

# To run GC-MS / LC-MS variants with the "wide" loader, point parent_dir /
# output_dir / experiment_name at the relevant assay, e.g.:
#   parent_dir = ".../proteomics/GC_MS", experiment_name = "GC_MS"
#   parent_dir = ".../proteomics/LC_MS_neg", experiment_name = "LC_MS_neg"
#   parent_dir = ".../proteomics/LC_MS_pos", experiment_name = "LC_MS_pos"

# =============================================================================
# 03_lc_feature_exploration.R  (one-off exploration, not a routine pipeline)
# =============================================================================
LC_EXPLORE <- list(
  output_dir = "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_pos/clustered_89",
  experiment_name = "LC_MS_pos_89",

  # Explicit input file(s) (sample x feature wide CSVs).
  files = c("/data/local/jy1008/SaMu/proteomics/LC_MS_pos/data/LCMS_processed_blocks8,9,1PY2_SaMu_298_371and5PY2_25PY2_pos.csv"),

  # Keep peaks present (non-zero) in at least this fraction of samples.
  peak_filter_threshold = 0.25,
  # Correlation threshold for building the peak co-occurrence network.
  corr_threshold = 0.8,
  # Drop correlation-network clusters smaller than this many peaks.
  min_cluster_size = 2,
  # Collapse each cluster to a representative value: "mean" or "max".
  cluster_summary = "mean"
)
