# =============================================================================
# config.R — central configuration for the metagenomics downstream analysis
# (MetaPhlAn profiles -> metadata join -> filtering -> diversity -> DESeq2).
#
# Source this at the top of every numbered script:
#   source("config.R")
#
# Edit the values here to point the analysis at a new study. The numbered
# scripts (01_*, 02_*, 03_*) should not need study-specific edits.
# =============================================================================

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# Directory of MetaPhlAn per-sample profiles (*_profile.txt).
INPUT_DIR <- "/data/local/jy1008/SaMu/metaphlan_out_02102026/all_merged_fastqs"

# Output / working directory. All results (CSVs, plots, .rds) are written here.
OUTPUT_DIR <- "/data/local/jy1008/SaMu/results/latest/metagenomics_R"

# Sample metadata CSV.
METADATA_CSV <- "/data/local/jy1008/SaMu/metadata/SaMu_sarcopeniestatus_majorcovariates_v13_16012026.csv"

# ---------------------------------------------------------------------------
# Output filenames
# ---------------------------------------------------------------------------
# A single tag appended to dated/versioned outputs so a new run does not
# overwrite a previous one. Set to "" for no suffix.
RUN_TAG <- "03132026"

# Intermediates passed between the numbered scripts (written into OUTPUT_DIR).
META_FILTERED_RDS <- "meta_filtered.rds"   # long-format filtered species table
META_DF_RDS       <- "meta_df.rds"         # per-sample metadata (aligned)

# ---------------------------------------------------------------------------
# Metadata schema
# ---------------------------------------------------------------------------
# Column in METADATA_CSV that uniquely identifies a subject/record; used as the
# metadata row name.
RECORD_ID_COL <- "record_id"

# Restrict to a cohort with a flag column == 1 (set FULLSAMU_COL <- NA to skip).
FULLSAMU_COL <- "Full.SaMu"

# Columns required to be complete for a sample to enter the models. A sample
# missing any of these is dropped.
COLS_OF_INTEREST <- c("sarc_status_bin", "age_def", "sex", "smke",
                      "nutr_score", "bmi", "alco")

# Columns to retain from the (wide) metadata before analysis. File_ID and
# sarc_status_bin are added by the recode step, so they are listed here too.
META_KEEP_COLS <- c(
  "record_id", "Block", "Full.SaMu", "sarc_status",
  "EWGSOP_strength", "EWGSOP_mass", "EWGSOP_performance",
  "age_def", "sex", "smke", "alco", "nutr_score", "bmi", "stvol",
  "File_ID", "sarc_status_bin"
)

# ---------------------------------------------------------------------------
# File_ID construction
# ---------------------------------------------------------------------------
# Each sequencing block stores its raw-file reference in a different metadata
# column, and the sample ID is embedded in that string with a block-specific
# format. For each block, give the column name and a regex that extracts the
# sample ID matching the MetaPhlAn output basename (<File_ID>_profile.txt).
#
# Processed top to bottom; a non-empty cell in a later block overrides an
# earlier one for the same row. Add or remove entries as blocks change.
FILE_ID_BLOCKS <- list(
  list(column = "block1.20231011_AV224503_4520_1.RawData.4520.tar",
       regex  = "SaMu[0-9]+"),
  list(column = "block2.20240327_AV224503_4734_1.RawData.4734.tar",
       regex  = "\\d+_\\d+"),
  list(column = "block5.20250124_AV242402_4915_1.RawData.4915.tar",
       regex  = "\\d+_\\d+_\\d+_SaMu\\d+"),
  list(column = "block89a.2025024_AV24242_5026_5092_1.RawData.5092",
       regex  = "\\d+_\\d+_Libr\\d+_SaMu\\d+")
)

# ---------------------------------------------------------------------------
# Abundance filtering
# ---------------------------------------------------------------------------
# A species is kept if its relative abundance (%) is >= ABUND_THRESHOLD in at
# least MIN_SAMPLES samples. Relative abundances are out of 100, so 0.1 = 0.1%.
ABUND_THRESHOLD <- 0.1
MIN_SAMPLES     <- 10

# Elbow plot sweep (diagnostic for choosing the two values above).
ELBOW_ABUND_THRESHOLDS    <- seq(0, 2, by = 0.02)
ELBOW_PREVALENCE_THRESHOLDS <- c(1, 5, 10, 20)

# Family/genus stacked-barplot: collapse taxa below this mean abundance (%)
# into "Other"; "unclassified" taxa are always collapsed.
BARPLOT_MEAN_ABUND_CUTOFF <- 1.5

# ---------------------------------------------------------------------------
# Grouping variable and differential-abundance model
# ---------------------------------------------------------------------------
# Binary grouping used throughout (diversity tests, ordination color, DESeq2).
GROUP_VAR     <- "sarc_status_bin"
GROUP_LEVELS  <- c("NoSarc", "Sarc")   # first level is the reference

# Continuous covariates to z-score before DESeq2 (a "<col>_scaled" column is
# created for each).
SCALE_COVARIATES <- c("age_def", "nutr_score", "bmi")

# DESeq2 design (right-hand side). Use the *_scaled names for scaled covariates
# and put the grouping variable last so its coefficient is the contrast.
DESEQ_DESIGN <- ~ age_def_scaled + sex + smke + alco +
                  nutr_score_scaled + bmi_scaled + sarc_status_bin

# Coefficient passed to lfcShrink (must match a resultsNames(dds) entry).
DESEQ_COEF <- "sarc_status_bin_Sarc_vs_NoSarc"

# Significance cutoff for the DESeq2 dotplot.
PADJ_CUTOFF <- 0.05

# Plots split species into two facets at this absolute log2 fold change.
LFC_FACET_CUTOFF <- 2
