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
  library(jsonlite)
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
# find_density_valley — locate the local minimum in a kernel density estimate
# of log10(x) sitting between its two tallest, sufficiently-separated peaks.
# Used to find a data-driven filtering threshold when a diagnostic histogram
# (e.g. abundance_histogram_mean.png) looks bimodal — the valley between two
# real populations of features is a more defensible cutoff than a threshold
# borrowed from a different tool/dataset (see abundance_histogram_mean.png
# discussion: MIN_MEAN_ABUNDANCE from MaAsLin2's tutorial can land anywhere
# relative to this dataset's actual distribution, including right on a peak).
#
# min_peak_separation is in log10 units — guards against picking two small
# adjacent bumps from density-estimation noise as if they were the real modes.
# Returns NA (with a warning) if fewer than two well-separated peaks are found.
# ---------------------------------------------------------------------------
find_density_valley <- function(x, min_peak_separation = 0.5) {
  x <- x[!is.na(x) & x > 0]
  log_x <- log10(x)
  d <- density(log_x, n = 1024)

  is_peak <- c(FALSE, diff(sign(diff(d$y))) == -2, FALSE)
  peak_idx <- which(is_peak)
  if (length(peak_idx) < 2) {
    warning("Fewer than two local maxima found in the density estimate — ",
            "distribution may not be bimodal.")
    return(NA_real_)
  }

  peaks <- data.frame(x = d$x[peak_idx], y = d$y[peak_idx]) %>%
    arrange(desc(y))
  peak1 <- peaks[1, ]
  candidates <- peaks[abs(peaks$x - peak1$x) >= min_peak_separation, ]
  if (nrow(candidates) == 0) {
    warning("No second peak found at least ", min_peak_separation,
            " log10 units from the tallest peak — distribution may not be cleanly bimodal.")
    return(NA_real_)
  }
  peak2 <- candidates[1, ]

  lo <- min(peak1$x, peak2$x); hi <- max(peak1$x, peak2$x)
  between <- d$x >= lo & d$x <= hi
  valley_x <- d$x[between][which.min(d$y[between])]
  10 ^ valley_x
}

# ---------------------------------------------------------------------------
# read_fastp_depth — read all *.json fastp reports in a directory, extract
# each sample's post-filtering read count (needed as a MaAsLin3 covariate —
# see config.R's FIXED_EFFECTS comment for why), and resolve each filename to
# a File_ID using the same extract_sample_id() logic already used for HUMAnN
# sample columns. fastp JSONs come from the same fastq-derived naming
# convention (3_remove_host_reads.sh / 1_combine_L001_L002.sh etc.), so the
# same STRIP_SUFFIX_REGEX / FILE_ID_REGEXES apply without modification.
#
# Returns a tibble: File_ID, read_depth.
#
# NOTE: fastp's JSON schema stores the post-QC read count at
# summary$after_filtering$total_reads (R1+R2 combined for paired-end, per
# fastp's documented schema). Verify this against one of your actual files
# before trusting it — e.g. str(jsonlite::fromJSON(files[1])$summary) — since
# schema details can vary slightly by fastp version/invocation flags.
# ---------------------------------------------------------------------------
read_fastp_depth <- function(dir) {
  files <- list.files(dir, pattern = "\\.json$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No .json files found in ", dir, " — check FASTP_JSON_DIR in config.R.")
  }
  ids <- extract_sample_id(file_path_sans_ext(basename(files)))
  depths <- vapply(files, function(f) {
    j <- fromJSON(f)
    as.numeric(j$summary$after_filtering$total_reads)
  }, numeric(1), USE.NAMES = FALSE)
  tibble(File_ID = ids, read_depth = depths)
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