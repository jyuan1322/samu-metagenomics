# =============================================================================
# 03_maaslin.R
#
# Run MaAsLin3 on the filtered HUMAnN feature table (from
# 01_load_and_filter.R) against FIXED_EFFECTS / RANDOM_EFFECTS (config.R).
#
# MaAsLin3 (not MaAsLin2): current biobakery-supported version, tests both
# abundance and prevalence associations, better handles compositionality.
#   BiocManager::install("biobakery/maaslin3")
#
# Run from this directory, after 01_load_and_filter.R:  Rscript 03_maaslin.R
# =============================================================================
source("config.R")
source("utils.R")

suppressPackageStartupMessages(library(maaslin3))

setwd(OUTPUT_DIR)

feature_mat  <- readRDS(FILTERED_FEATURES_RDS)
meta_aligned <- readRDS(META_ALIGNED_RDS)

message(sprintf("Loaded filtered feature table: %d features x %d samples",
                nrow(feature_mat), ncol(feature_mat)))

# read depth is added in 01_load_and_filter.R as a column in meta_aligned
# read_depth <- read_depth_df$read_depth[match(meta_aligned$File_ID, read_depth_df$File_ID)]
# meta_aligned$read_depth <- read_depth

# ---------------------------------------------------------------------------
# Metadata for MaAsLin3: rownames = File_ID, columns = model variables.
# sarc_status_bin arrives already recoded + factored (levels NoSarc, Sarc)
# from metagenomics_R's meta_df.rds — the loop below only coerces as a
# fallback (e.g. for a future numeric/low-cardinality fixed effect added to
# FIXED_EFFECTS that hasn't been explicitly factored upstream).
# ---------------------------------------------------------------------------
model_vars <- c(FIXED_EFFECTS, RANDOM_EFFECTS)
meta_maaslin <- as.data.frame(meta_aligned)
rownames(meta_maaslin) <- meta_maaslin$File_ID
meta_maaslin <- meta_maaslin[, model_vars, drop = FALSE]

for (fx in FIXED_EFFECTS) {
  if (!is.factor(meta_maaslin[[fx]]) &&
      length(unique(na.omit(meta_maaslin[[fx]]))) <= 5) {
    meta_maaslin[[fx]] <- as.factor(meta_maaslin[[fx]])
    message(sprintf("Coerced '%s' to factor (levels: %s) — expected this to already be a factor; double check upstream recoding if so.",
                    fx, paste(levels(meta_maaslin[[fx]]), collapse = ", ")))
  }
}

# feature_mat is features x samples; MaAsLin3 wants samples as rows to line
# up 1:1 with meta_maaslin.
feature_df <- as.data.frame(t(feature_mat))
stopifnot(identical(rownames(feature_df), rownames(meta_maaslin)))

# ---------------------------------------------------------------------------
# Run MaAsLin3
# ---------------------------------------------------------------------------
out_dir <- tag_filename(paste0("maaslin3_", FEATURE_TABLE))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf("Running MaAsLin3 (fixed_effects = %s%s) -> %s",
                paste(FIXED_EFFECTS, collapse = ", "),
                if (length(RANDOM_EFFECTS)) paste0(", random_effects = ", paste(RANDOM_EFFECTS, collapse = ", ")) else "",
                out_dir))

fit_data <- maaslin3(
  input_data        = feature_df,
  input_metadata    = meta_maaslin,
  output            = out_dir,
  fixed_effects     = FIXED_EFFECTS,
  random_effects    = if (length(RANDOM_EFFECTS)) RANDOM_EFFECTS else NULL,
  normalization     = MAASLIN_NORMALIZATION,
  transform         = MAASLIN_TRANSFORM,
  max_significance  = MAASLIN_MAX_SIGNIFICANCE,
  augment           = TRUE,
  standardize       = TRUE,
  plot_summary_plot = TRUE,
  plot_associations = TRUE
)

message(sprintf("03_maaslin.R complete. MaAsLin3 output in: %s", file.path(OUTPUT_DIR, out_dir)))
message("Key files: significant_results.tsv (FDR-significant abundance + prevalence associations),")
message("           all_results.tsv (every feature x fixed-effect term tested)")
