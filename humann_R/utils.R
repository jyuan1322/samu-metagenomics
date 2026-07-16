# =============================================================================
# utils.R — helper functions for the HUMAnN downstream analysis.
# Source after config.R:  source("utils.R")
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tools)
  library(tibble)
})

# ---------------------------------------------------------------------------
# tag_filename — identical convention to metagenomics_R/utils.R: insert
# RUN_TAG before the extension of an output filename.
#   tag_filename("maaslin_input.csv") -> "maaslin_input_07162026.csv"
# ---------------------------------------------------------------------------
tag_filename <- function(name) {
  if (is.null(RUN_TAG) || RUN_TAG == "") return(name)
  ext  <- file_ext(name)
  stem <- file_path_sans_ext(name)
  if (ext == "") paste0(stem, "_", RUN_TAG)
  else paste0(stem, "_", RUN_TAG, ".", ext)
}

# ---------------------------------------------------------------------------
# read_humann_table — reads a humann_join_tables + humann_renorm_table
# output (features x samples, tab-delimited). First column is the feature ID
# ("# Pathway" or "# Gene Family"); every other column is one sample.
#
# Uses data.table::fread rather than readr::read_tsv here specifically
# because joined_genefamilies*.tsv runs ~3GB (per 7_merge_humann.sh output) —
# fread is dramatically faster/lower-memory for a file this size. Returns a
# tibble so everything downstream can stay in the usual dplyr style.
# ---------------------------------------------------------------------------
read_humann_table <- function(path) {
  df <- data.table::fread(path, sep = "\t", header = TRUE, check.names = FALSE) %>%
    as_tibble()
  names(df)[1] <- "feature_id"
  df
}

# ---------------------------------------------------------------------------
# extract_sample_id — resolve a HUMAnN sample column name down to the same
# File_ID format used in metadata.
#
#   1. Strip known unit suffixes (STRIP_SUFFIX_REGEX): "..._Abundance" etc.
#   2. Try each FILE_ID_REGEXES pattern in turn (same block-specific patterns
#      as FILE_ID_BLOCKS in metagenomics_R/config.R, since File_ID is
#      embedded the same way in both MetaPhlAn's and HUMAnN's output names,
#      just with extra pipeline-stage suffixes on the HUMAnN side — see the
#      comment above FILE_ID_REGEXES in this pipeline's config.R). First
#      pattern to match wins, so order matters where one pattern is a prefix
#      of another (block2's "\d+_\d+" is a prefix of block5/89a's patterns —
#      keep it listed last, same ordering caveat noted for FILE_ID_BLOCKS in
#      metagenomics_R/config.R).
#   3. If nothing matches, fall back to the stripped string as-is and rely on
#      diagnose_sample_matching() to surface the miss.
# ---------------------------------------------------------------------------
extract_sample_id <- function(names_vec) {
  stripped <- str_remove(names_vec, STRIP_SUFFIX_REGEX)
  extracted <- rep(NA_character_, length(stripped))
  for (pattern in FILE_ID_REGEXES) {
    still_unmatched <- is.na(extracted)
    if (!any(still_unmatched)) break
    hit <- str_extract(stripped[still_unmatched], pattern)
    extracted[still_unmatched] <- hit
  }
  ifelse(is.na(extracted), stripped, extracted)
}

# ---------------------------------------------------------------------------
# diagnose_sample_matching — prints match counts and a few examples from each
# side that failed to match, to make it fast to spot a wrong regex instead of
# silently dropping samples. Call this before trusting alignment counts.
# ---------------------------------------------------------------------------
diagnose_sample_matching <- function(humann_sample_ids, file_ids) {
  matched <- intersect(humann_sample_ids, file_ids)
  unmatched_humann <- setdiff(humann_sample_ids, file_ids)
  unmatched_meta   <- setdiff(file_ids, humann_sample_ids)

  message(sprintf(
    "Sample matching: %d matched, %d HUMAnN samples unmatched, %d metadata File_IDs unmatched",
    length(matched), length(unmatched_humann), length(unmatched_meta)
  ))
  if (length(unmatched_humann) > 0) {
    message("  Example unmatched HUMAnN sample IDs (after extraction): ",
            paste(head(unmatched_humann, 5), collapse = ", "))
  }
  if (length(unmatched_meta) > 0) {
    message("  Example unmatched metadata File_IDs: ",
            paste(head(unmatched_meta, 5), collapse = ", "))
  }
  invisible(list(matched = matched, unmatched_humann = unmatched_humann,
                 unmatched_meta = unmatched_meta))
}

# ---------------------------------------------------------------------------
# split_unstratified — HUMAnN feature IDs use "|" to mark species-stratified
# rows, e.g. "PWY-101|g__Bacteroides.s__fragilis". Community-total rows have
# no "|". Keep only unstratified rows for a per-pathway/per-gene-family
# association test.
# ---------------------------------------------------------------------------
split_unstratified <- function(df) {
  df %>% filter(!str_detect(feature_id, "\\|"))
}

# ---------------------------------------------------------------------------
# feature_stats — per-feature mean relative abundance, prevalence, and
# coefficient of variation, all computed over nonzero values only (so CV
# reflects genuine variability rather than being dominated by how often a
# feature is simply absent).
# mat: numeric matrix, features x samples, rownames = feature_id.
# ---------------------------------------------------------------------------
feature_stats <- function(mat) {
  n_samples <- ncol(mat)
  stats <- apply(mat, 1, function(x) {
    # apply(MARGIN=1) passes each row in with names(x) set to colnames(mat)
    # (the sample IDs). Strip them here — otherwise, when a feature has
    # exactly one nonzero sample, nz below is a length-1 *named* vector, and
    # c(mean_abund = nz) doesn't produce a vector named "mean_abund" as
    # intended — R concatenates the names into "mean_abund.<samplename>"
    # instead, since nz already carries a name. That breaks column naming
    # once results from all features get combined into one matrix/data frame.
    x <- unname(x)
    nz <- x[x > 0]
    prevalence <- length(nz) / n_samples
    if (length(nz) < 2) {
      return(c(mean_abund = if (length(nz) == 1) nz else NA_real_,
               cv = NA_real_, prevalence = prevalence))
    }
    m <- mean(nz)
    s <- sd(nz)
    c(mean_abund = m, cv = if (m > 0) s / m else NA_real_, prevalence = prevalence)
  })
  as.data.frame(t(stats)) %>%
    rownames_to_column("feature_id")
}

# ---------------------------------------------------------------------------
# filter_by_abundance_cv — applies MIN_MEAN_ABUNDANCE, MIN_PREVALENCE, and
# MIN_CV together. Returns the filtered matrix plus the full stats table.
# ---------------------------------------------------------------------------
filter_by_abundance_cv <- function(mat, min_abund, min_prev, min_cv) {
  st <- feature_stats(mat)
  keep <- st$feature_id[
    !is.na(st$mean_abund) & st$mean_abund >= min_abund &
    st$prevalence >= min_prev &
    !is.na(st$cv) & st$cv >= min_cv
  ]
  list(
    mat = mat[rownames(mat) %in% keep, , drop = FALSE],
    stats = st,
    kept_ids = keep
  )
}

# ---------------------------------------------------------------------------
# df_to_matrix — feature_id column -> rownames, rest to numeric matrix.
# ---------------------------------------------------------------------------
df_to_matrix <- function(df) {
  m <- as.matrix(df %>% select(-feature_id))
  rownames(m) <- df$feature_id
  storage.mode(m) <- "double"
  m
}