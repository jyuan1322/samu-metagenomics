# =============================================================================
# utils.R — helper functions shared across the analysis scripts.
# Source after config.R:  source("utils.R")
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tools)
})

# ---------------------------------------------------------------------------
# build_file_id — populate a File_ID column on the metadata frame by extracting
# a sample ID from a block-specific column with a block-specific regex.
#
# Generalizes the original hardcoded four-block sequence: it loops over the
# FILE_ID_BLOCKS list from config.R. Blocks are applied in order, so a
# non-empty match in a later block overrides an earlier one for the same row.
#
# Args:
#   df     : metadata data.frame (wide, one row per record)
#   blocks : list of list(column=, regex=) entries (default: FILE_ID_BLOCKS)
# Returns: df with a populated character File_ID column.
# ---------------------------------------------------------------------------
build_file_id <- function(df, blocks = FILE_ID_BLOCKS) {
  df$File_ID <- NA_character_
  for (blk in blocks) {
    col <- blk$column
    if (!col %in% names(df)) {
      warning(sprintf("File_ID block column not found, skipping: %s", col))
      next
    }
    nonempty_idx <- trimws(df[[col]]) != "" & !is.na(df[[col]])
    df$File_ID[nonempty_idx] <- str_extract(df[[col]][nonempty_idx], blk$regex)
  }
  df
}

# ---------------------------------------------------------------------------
# recode_metadata — apply the study's covariate recoding. Edit here if the
# metadata coding scheme changes. Kept as one documented function so the
# encoding assumptions live in a single place.
# ---------------------------------------------------------------------------
recode_metadata <- function(df) {
  # smoking: any of {1,2,3} -> ever-smoker (1); 0 -> never (0); else NA
  df$smke <- ifelse(df$smke %in% c("1", "2", "3"), 1,
                    ifelse(df$smke == "0", 0, NA))
  # alcohol: {1,2} -> 0 ; {3,4} -> 1 ; else NA
  df$alco <- ifelse(df$alco %in% c("1", "2"), 0,
                    ifelse(df$alco %in% c("3", "4"), 1, NA))
  # sex: 0 -> M, 1 -> F, else NA
  df$sex <- ifelse(df$sex == "0", "M",
                   ifelse(df$sex == "1", "F", NA))

  df$age_def    <- as.numeric(df$age_def)
  df$bmi        <- as.numeric(df$bmi)
  df$nutr_score <- as.numeric(df$nutr_score)

  message("Converting sarc_status to numeric; non-numeric values become NA")
  df$sarc_status <- suppressWarnings(as.numeric(df$sarc_status))

  # Binarize the grouping variable: > 0 -> Sarc, else NoSarc.
  df$sarc_status_bin <- ifelse(df$sarc_status > 0, "Sarc", "NoSarc")
  df$sarc_status_bin <- factor(df$sarc_status_bin, levels = GROUP_LEVELS)
  df
}

# ---------------------------------------------------------------------------
# read_metaphlan — read one MetaPhlAn profile, returning a tidy long frame
# with Sample, clade_name, relative_abundance, and estimated read counts.
# ---------------------------------------------------------------------------
read_metaphlan <- function(file) {
  lines <- readLines(file)
  header_line <- lines[grepl("^#clade_name", lines)]
  header <- strsplit(sub("^#", "", header_line), "\\s+")[[1]]

  df <- read_tsv(file, comment = "#", col_names = header, col_types = cols())
  df$Sample <- file_path_sans_ext(basename(file))

  df %>%
    select(Sample, clade_name, relative_abundance,
           estimated_number_of_reads_from_the_clade)
}

# ---------------------------------------------------------------------------
# tag_filename — insert RUN_TAG before the extension of an output filename.
#   tag_filename("deseq2_results.csv") -> "deseq2_results_03132026.csv"
# Returns the name unchanged when RUN_TAG is "".
# ---------------------------------------------------------------------------
tag_filename <- function(name) {
  if (is.null(RUN_TAG) || RUN_TAG == "") return(name)
  ext  <- file_ext(name)
  stem <- file_path_sans_ext(name)
  if (ext == "") paste0(stem, "_", RUN_TAG)
  else paste0(stem, "_", RUN_TAG, ".", ext)
}

# ---------------------------------------------------------------------------
# renormalize_relab — rescale each sample's relative abundances to sum to 100.
# ---------------------------------------------------------------------------
renormalize_relab <- function(df, sample_col = "Sample",
                              value_col = "relative_abundance") {
  df %>%
    group_by(.data[[sample_col]]) %>%
    mutate(!!value_col := .data[[value_col]] /
             sum(.data[[value_col]]) * 100) %>%
    ungroup()
}

# Consistent palette for the stacked barplots (10 distinct + grey for "Other").
BARPLOT_COLORS <- c(
  "#0074D9", "#FF851B", "#2ECC40", "#FF4136", "#B10DC9",
  "#FFDC00", "#39CCCC", "#85144B", "#01FF70", "#001F3F", "#808080"
)
