# =============================================================================
# 04_ordinal_lasso.R
#
# Standalone L1-penalized (lasso) ordinal regression on the ORIGINAL ordered
# sarc_status outcome (0 = control, 1/2/3 = increasing sarcopenia severity),
# using a parallel cumulative logit (proportional-odds) model fit by ordinalNet.
#
# This is a separate analysis track from the binary pipeline: the other models
# (lasso/RF/NN) predict the binarized sarc_status_bin, whereas this one keeps
# all ordered levels of sarc_status. It therefore reads the metadata fresh and
# recodes with extreme_cases_only = FALSE, so it is unaffected by the
# EXTREME_CASES_ONLY setting used elsewhere (which would otherwise drop the
# sarc_status == 1 level this model needs).
#
# Method notes (per Wurm, Hanlon & Rathouz 2021, the ordinalNet paper):
#   * family = "cumulative", link = "logit", parallelTerms = TRUE,
#     nonparallelTerms = FALSE  -> proportional-odds model.
#   * alpha = 1                 -> pure lasso penalty.
#   * predictors are scaled to unit variance for balanced penalization.
#   * ordinalNet uses an inverse sign convention; fitted slope coefficients are
#     multiplied by -1 here to match the conventional direction (positive =>
#     higher feature value associates with higher severity).
#   * lambda is tuned by K-fold CV on out-of-sample log-likelihood
#     (ordinalNetTune); ordinalNetCV gives a cross-validated performance
#     estimate of the tuned model.
#
# Run:  Rscript 04_ordinal_lasso.R
# =============================================================================
source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ordinalNet)
})

setwd(OUTPUT_DIR)

set.seed(ORDINAL_SEED)

# ---------------------------------------------------------------------------
# Features: CLR-transformed species relative abundances (Sample x Species),
# read from the filtered long table written by 01_load_and_filter.R.
# ---------------------------------------------------------------------------
meta_filtered <- readRDS(META_FILTERED_RDS)
meta_filtered$Sample <- gsub("_profile", "", meta_filtered$Sample)

relab_wide <- meta_filtered %>%
  select(Sample, Species, relative_abundance) %>%
  pivot_wider(names_from = Species, values_from = relative_abundance,
              values_fill = 0)

sample_ids <- relab_wide$Sample
feat_mat   <- as.matrix(relab_wide[, -1])
rownames(feat_mat) <- sample_ids

# CLR transform: add a pseudo-count, renormalize, log, subtract the per-sample
# geometric mean (in log space). Matches the metagenomics Python loader.
clr_transform <- function(x, pseudo_count = ORDINAL_PSEUDO_COUNT) {
  x <- x + pseudo_count
  x <- x / rowSums(x)
  log_x <- log(x)
  log_x - rowMeans(log_x)
}
clr_mat <- clr_transform(feat_mat)

# ---------------------------------------------------------------------------
# Outcome: the ORIGINAL ordered sarc_status (0/1/2/3), recoded fresh from the
# metadata CSV with extreme_cases_only = FALSE so every level is retained.
# ---------------------------------------------------------------------------
meta_raw <- read.csv(METADATA_CSV, header = TRUE, stringsAsFactors = FALSE)
meta_raw <- build_file_id(meta_raw)

if (!is.na(FULLSAMU_COL)) {
  meta_raw <- meta_raw[meta_raw[[FULLSAMU_COL]] == 1, ]
}

meta_raw <- recode_metadata(meta_raw, extreme_cases_only = FALSE)

# Drop subjects with no sequencing File_ID (present for other assays only).
meta_raw <- meta_raw[!is.na(meta_raw$File_ID), ]
if (anyDuplicated(meta_raw$File_ID)) {
  dups <- unique(meta_raw$File_ID[duplicated(meta_raw$File_ID)])
  stop("Duplicate File_ID(s) after build_file_id: ",
       paste(dups, collapse = ", "),
       " — check FILE_ID_BLOCKS regexes in config.R")
}
rownames(meta_raw) <- meta_raw$File_ID

# ---------------------------------------------------------------------------
# Align features and metadata; scale covariates and features to unit variance.
# ---------------------------------------------------------------------------
common <- intersect(rownames(clr_mat), rownames(meta_raw))
clr_mat  <- clr_mat[common, , drop = FALSE]
meta_sub <- meta_raw[common, ]

# Covariate block (continuous covariates kept numeric; sex -> 0/1; binary
# smke/alco already 0/1 from recode_metadata).
meta_sub$sex_bin <- ifelse(meta_sub$sex == "F", 1,
                           ifelse(meta_sub$sex == "M", 0, NA))
cov_cols <- c(ORDINAL_CONT_COVARIATES, "sex_bin", ORDINAL_BIN_COVARIATES)
cov_mat  <- as.matrix(meta_sub[, cov_cols])

# Outcome as an ordered factor over the full 0/1/2/3 range.
y_raw <- meta_sub$sarc_status

# Optionally pool sparse top levels into a single ceiling level (e.g. 2+).
collapse_label <- NA_character_
if (!is.na(ORDINAL_MAX_LEVEL)) {
  collapse_label <- paste0(ORDINAL_MAX_LEVEL, "+")
  y_raw <- ifelse(!is.na(y_raw) & y_raw >= ORDINAL_MAX_LEVEL,
                  ORDINAL_MAX_LEVEL, y_raw)
}

# Complete-case filter over features, covariates, and outcome.
keep <- stats::complete.cases(clr_mat) &
        stats::complete.cases(cov_mat) &
        !is.na(y_raw)
clr_mat <- clr_mat[keep, , drop = FALSE]
cov_mat <- cov_mat[keep, , drop = FALSE]
y_raw   <- y_raw[keep]

x_mat <- cbind(cov_mat, clr_mat)

# Scale all predictors to unit variance (drop zero-variance columns, which
# ordinalNet cannot scale and which carry no signal).
nzv <- apply(x_mat, 2, function(col) stats::sd(col, na.rm = TRUE) > 0)
if (any(!nzv)) {
  message(sum(!nzv), " zero-variance predictor(s) dropped: ",
          paste(colnames(x_mat)[!nzv], collapse = ", "))
}
x_mat <- x_mat[, nzv, drop = FALSE]
x_scaled <- scale(x_mat)

# Build the ordered factor. When collapsing, the top numeric level is
# relabelled (e.g. "2+") so the pooling is explicit in the outputs.
ord_levels <- sort(unique(y_raw))
ord_labels <- as.character(ord_levels)
if (!is.na(ORDINAL_MAX_LEVEL)) {
  ord_labels[ord_levels == ORDINAL_MAX_LEVEL] <- collapse_label
  message("Collapsed sarc_status >= ", ORDINAL_MAX_LEVEL,
          " into level '", collapse_label, "'")
}
y_ord <- factor(y_raw, levels = ord_levels, labels = ord_labels, ordered = TRUE)
message("Ordinal outcome level counts (sarc_status):")
print(table(y_ord))
message(nrow(x_scaled), " samples x ", ncol(x_scaled), " predictors")

# ---------------------------------------------------------------------------
# Tune lambda by K-fold CV (best out-of-sample log-likelihood), then refit.
# ---------------------------------------------------------------------------
tune_fit <- ordinalNetTune(
  x_scaled, y_ord,
  family = "cumulative", link = "logit",
  parallelTerms = TRUE, nonparallelTerms = FALSE,
  alpha = 1,
  nFolds = ORDINAL_NFOLDS,
  printProgress = FALSE
)

# Best lambda = highest mean out-of-sample log-likelihood across folds.
best_idx <- which.max(rowMeans(tune_fit$loglik))
message("Selected lambda index ", best_idx, " of ",
        length(tune_fit$fit$lambdaVals))

# ---------------------------------------------------------------------------
# Cross-validated performance estimate of the tuned model.
# ---------------------------------------------------------------------------
cv_fit <- ordinalNetCV(
  x_scaled, y_ord,
  family = "cumulative", link = "logit",
  parallelTerms = TRUE, nonparallelTerms = FALSE,
  alpha = 1,
  nFolds = ORDINAL_NFOLDS,
  tuneMethod = "cvLoglik",
  printProgress = FALSE
)
cv_summary <- summary(cv_fit)
write.csv(cv_summary, tag_filename("ordinal_lasso_cv_performance.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Coefficients at the selected lambda. ordinalNet uses an inverse sign
# convention, so slope coefficients are multiplied by -1. The parallel model
# has one shared slope vector plus (K-1) intercepts.
# ---------------------------------------------------------------------------
coef_mat <- coef(tune_fit$fit, whichLambda = best_idx, matrix = TRUE)

# Intercept rows are the cumulative-logit thresholds; the remaining rows are
# the shared predictor slopes. Identify slopes by predictor name.
all_terms <- rownames(coef_mat)
slope_rows <- all_terms %in% colnames(x_scaled)

slopes <- coef_mat[slope_rows, 1]
slopes <- -slopes  # invert sign convention -> conventional direction

coef_df <- data.frame(
  feature = names(slopes),
  coef = as.numeric(slopes),
  odds_ratio = exp(as.numeric(slopes)),
  selected = as.numeric(slopes) != 0,
  row.names = NULL
)
coef_df <- coef_df[order(-abs(coef_df$coef)), ]

write.csv(coef_df, tag_filename("ordinal_lasso_coefficients.csv"),
          row.names = FALSE)

# Thresholds (intercepts) saved separately for completeness.
intercept_df <- data.frame(
  threshold = all_terms[!slope_rows],
  value = as.numeric(coef_mat[!slope_rows, 1]),
  row.names = NULL
)
write.csv(intercept_df, tag_filename("ordinal_lasso_thresholds.csv"),
          row.names = FALSE)

n_selected <- sum(coef_df$selected)
message(n_selected, " of ", nrow(coef_df),
        " predictors retained (non-zero) at the selected lambda")
message("Wrote: ",
        tag_filename("ordinal_lasso_coefficients.csv"), ", ",
        tag_filename("ordinal_lasso_cv_performance.csv"), ", ",
        tag_filename("ordinal_lasso_thresholds.csv"))