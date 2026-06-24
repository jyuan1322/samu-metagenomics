# =============================================================================
# metadata_utils.R
#
# Metadata loading, covariate recoding, and cohort filtering for the
# other-omics analyses. These scripts recode from the raw metadata CSV rather
# than reusing the metagenomics pipeline's output, because the other-omics
# assays can include samples that were dropped from (or never entered) the
# metagenomics cohort. Keeping the recode in one place here avoids it drifting
# between the omics scripts (e.g. one binarizing smoking/alcohol while another
# assumes it is already binarized).
#
# Source after the caller's config has defined GROUP_LEVELS (and, for
# load_samu_metadata, the relevant column names). Functions take their inputs
# as arguments rather than reading globals, so they are safe to source anywhere.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ---------------------------------------------------------------------------
# recode_samu_covariates — apply the standard SaMu covariate recoding.
#
# Args:
#   df            : metadata data.frame
#   group_levels  : c(reference, other) for the binarized grouping factor
#   recode_smoke_alcohol : TRUE to recode raw smke/alco codes; set FALSE when
#                          the metadata file already stores them binarized.
# Returns df with recoded smke/alco/sex, numeric sarc_status, and the
# sarc_status_bin factor.
#
# Coding scheme (when recode_smoke_alcohol = TRUE):
#   smke: {1,2,3} -> 1 (ever), 0 -> 0 (never), else NA
#   alco: {1,2}  -> 0,         {3,4} -> 1,     else NA
#   sex : 0 -> "M", 1 -> "F", else NA
#   sarc_status_bin: sarc_status > 0 -> "Sarc", else "NoSarc"
# ---------------------------------------------------------------------------
recode_samu_covariates <- function(df,
                                   group_levels = c("NoSarc", "Sarc"),
                                   recode_smoke_alcohol = TRUE) {
  if (recode_smoke_alcohol) {
    df$smke <- ifelse(df$smke %in% c("1", "2", "3"), 1,
                      ifelse(df$smke == "0", 0, NA))
    df$alco <- ifelse(df$alco %in% c("1", "2"), 0,
                      ifelse(df$alco %in% c("3", "4"), 1, NA))
  }
  if ("sex" %in% names(df)) {
    df$sex <- ifelse(df$sex == "0", "M",
                     ifelse(df$sex == "1", "F", df$sex))
  }

  message("Converting sarc_status to numeric; non-numeric values become NA")
  df$sarc_status <- suppressWarnings(as.numeric(df$sarc_status))

  df$sarc_status_bin <- ifelse(df$sarc_status > 0, "Sarc", "NoSarc")
  df$sarc_status_bin <- factor(df$sarc_status_bin, levels = group_levels)
  df
}

# ---------------------------------------------------------------------------
# load_samu_metadata — read the metadata CSV, restrict to the cohort flag,
# keep the requested columns, and recode covariates.
#
# Args:
#   path            : metadata CSV path
#   keep_cols       : columns to retain (NULL keeps all)
#   fullsamu_col    : cohort flag column kept == 1 (NA to skip)
#   group_levels    : passed through to recode_samu_covariates
#   recode_smoke_alcohol : passed through to recode_samu_covariates
#   add_label       : if TRUE, add label = paste0("SaMu", record_id)
# ---------------------------------------------------------------------------
load_samu_metadata <- function(path,
                               keep_cols = NULL,
                               fullsamu_col = "Full.SaMu",
                               group_levels = c("NoSarc", "Sarc"),
                               recode_smoke_alcohol = TRUE,
                               add_label = FALSE) {
  df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)

  if (add_label) df$label <- paste0("SaMu", df$record_id)

  if (!is.null(keep_cols)) {
    missing <- setdiff(keep_cols, names(df))
    if (length(missing) > 0)
      warning("Requested keep_cols not in metadata: ",
              paste(missing, collapse = ", "))
    df <- df[, intersect(keep_cols, names(df)), drop = FALSE]
  }

  if (!is.na(fullsamu_col) && fullsamu_col %in% names(df)) {
    df <- df %>% filter(.data[[fullsamu_col]] == 1)
  }

  recode_samu_covariates(df,
                         group_levels = group_levels,
                         recode_smoke_alcohol = recode_smoke_alcohol)
}
