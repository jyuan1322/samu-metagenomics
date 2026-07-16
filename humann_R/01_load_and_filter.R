# =============================================================================
# 01_load_and_filter.R
#
# Load the joined/renormalized HUMAnN table (from 7_merge_humann.sh), keep
# unstratified rows, align sample columns to File_ID using the already
# recoded/filtered metadata from metagenomics_R (METAGENOMICS_META_DF_RDS —
# no File_ID construction or covariate recoding happens here, that logic
# lives once in metagenomics_R/utils.R), and filter features by mean relative
# abundance + prevalence + coefficient of variation. Saves the intermediates
# consumed by 02_maaslin.R.
#
# Run from this directory:  Rscript 01_load_and_filter.R
# =============================================================================
source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(purrr); library(tidyr)
})

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
setwd(OUTPUT_DIR)

# ---------------------------------------------------------------------------
# Metadata: reuse the already-recoded, complete-case-filtered meta_df from
# metagenomics_R (File_ID + sarc_status_bin already built/recoded there).
# ---------------------------------------------------------------------------
meta_df <- readRDS(METAGENOMICS_META_DF_RDS)
message(sprintf("Loaded meta_df from metagenomics_R: %d subjects", nrow(meta_df)))
stopifnot("File_ID" %in% names(meta_df))
stopifnot(all(FIXED_EFFECTS %in% names(meta_df)))

# ---------------------------------------------------------------------------
# Feature table: load, keep unstratified rows, extract File_ID from sample
# column names
# ---------------------------------------------------------------------------
feature_path <- if (FEATURE_TABLE == "pathabundance") PATHABUNDANCE_FILE else GENEFAMILIES_FILE
message(sprintf("Loading %s: %s", FEATURE_TABLE, feature_path))

feat_df <- read_humann_table(feature_path)
feat_df <- split_unstratified(feat_df)
message(sprintf("%d unstratified features x %d samples (before filtering)",
                nrow(feat_df), ncol(feat_df) - 1))

mat <- df_to_matrix(feat_df)
colnames(mat) <- extract_sample_id(colnames(mat))

# ---------------------------------------------------------------------------
# Align to metadata by File_ID
# ---------------------------------------------------------------------------
diagnose_sample_matching(colnames(mat), meta_df$File_ID)

shared_ids <- intersect(colnames(mat), meta_df$File_ID)
if (length(shared_ids) == 0) {
  stop("No overlapping sample IDs between the HUMAnN table and meta_df$File_ID. ",
       "Check STRIP_SUFFIX_REGEX / FILE_ID_REGEXES in config.R against ",
       "colnames(read_humann_table(PATHABUNDANCE_FILE)) — see utils.R::extract_sample_id.")
}

mat <- mat[, shared_ids, drop = FALSE]
meta_aligned <- meta_df[match(shared_ids, meta_df$File_ID), ]
stopifnot(identical(colnames(mat), meta_aligned$File_ID))

# ---------------------------------------------------------------------------
# Feature stats (computed here, ahead of the elbow plots below, since the
# CV-filtered abundance elbow needs st$cv to know which features pass a given
# CV percentile).
# ---------------------------------------------------------------------------
st <- feature_stats(mat)
valid_cv <- st$cv[!is.na(st$cv)]

# MIN_CV is now derived from the data as a percentile (MIN_CV_PERCENTILE in
# config.R), not a hardcoded raw value — "keep the top P% most variable
# features," matching your advisor's convention. unname() strips the
# quantile()-added percentile label (e.g. "50%") so MIN_CV is a plain number
# for use in sprintf()/comparisons downstream.
MIN_CV <- unname(quantile(valid_cv, probs = 1 - MIN_CV_PERCENTILE / 100, na.rm = TRUE))
message(sprintf("MIN_CV derived from data: top %d%% -> CV >= %.4f", MIN_CV_PERCENTILE, MIN_CV))

# ---------------------------------------------------------------------------
# Elbow plot 2 (computed before plot 1 so plot 1 can reuse its CV cutoff):
# features retained vs. CV threshold, expressed as a percentile ("keep the
# top P% most variable features") rather than a raw CV value — matches your
# advisor's prior convention (e.g. "top 50%") and is dataset-relative, unlike
# a raw CV number tied to this dataset's particular scale/shape. One line per
# minimum-prevalence requirement.
# ---------------------------------------------------------------------------
elbow_cv_df <- expand_grid(
  top_pct = ELBOW_CV_PERCENTILES,
  prevalence_frac = ELBOW_PREVALENCE_FRACTIONS
) %>%
  mutate(FeaturesRetained = map2_int(top_pct, prevalence_frac, ~ {
    cv_cutoff <- quantile(valid_cv, probs = 1 - .x / 100, na.rm = TRUE)
    sum(!is.na(st$cv) & st$cv >= cv_cutoff & st$prevalence >= .y)
  }))

p_elbow_cv <- ggplot(elbow_cv_df,
                     aes(top_pct, FeaturesRetained, color = factor(prevalence_frac))) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = MIN_CV_PERCENTILE, color = "black", linetype = "dashed") +
  theme_minimal() +
  labs(x = "CV threshold, as \"keep top P% most variable\" (%)",
       y = sprintf("Number of %s features passing filter", FEATURE_TABLE),
       color = "Min. prevalence\n(fraction of samples)",
       title = "Elbow plot for CV filtering (percentile-based)")
ggsave(tag_filename("elbow_plot_cv.pdf"), p_elbow_cv, width = 10, height = 6)
# ggsave(tag_filename("elbow_plot_cv.png"), p_elbow_cv, width = 10, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Elbow plot 1: features retained vs. mean-abundance threshold, one line per
# minimum-prevalence requirement, drawn twice — once on the full feature set
# (solid lines) and once restricted to features passing the top
# MIN_CV_PERCENTILE% CV cutoff (dotted lines) — so you can see directly how
# much the CV filter shifts the abundance elbow, rather than inspecting the
# two filters in isolation.
#
# IMPORTANT: this must use the exact same per-feature quantities as
# filter_by_abundance_cv() — st$mean_abund (mean over NONZERO samples only)
# and st$prevalence (fraction nonzero, independent of any abundance
# threshold) — rather than a raw per-sample "value >= threshold" count.
# An earlier version of this plot swept rowSums(mat >= threshold) instead,
# which answers a different question ("in how many samples does the RAW
# value clear this threshold, jointly with a prevalence requirement") than
# what the actual filter checks ("is this feature's mean-of-nonzero-values
# above the threshold, evaluated independently of prevalence"). Those give
# different feature counts at the "same" threshold/prevalence_frac pair, so
# the plot wasn't a faithful preview of what filter_by_abundance_cv() would
# actually keep. Using st$mean_abund/st$prevalence directly (same as the CV
# elbow already does with st$cv) fixes that — and is also much cheaper, since
# it reuses the already-computed per-feature stats instead of re-scanning the
# full matrix once per threshold.
# ---------------------------------------------------------------------------
label_no_cv <- "No CV filter"
label_cv    <- sprintf("CV filter: top %d%%", MIN_CV_PERCENTILE)

elbow_abund_df <- expand_grid(
  threshold = ELBOW_ABUND_THRESHOLDS,
  prevalence_frac = ELBOW_PREVALENCE_FRACTIONS
) %>%
  mutate(FeaturesRetained = map2_int(threshold, prevalence_frac, ~ {
    sum(!is.na(st$mean_abund) & st$mean_abund >= .x & st$prevalence >= .y)
  }), cv_filter = label_no_cv)

elbow_abund_cv_df <- expand_grid(
  threshold = ELBOW_ABUND_THRESHOLDS,
  prevalence_frac = ELBOW_PREVALENCE_FRACTIONS
) %>%
  mutate(FeaturesRetained = map2_int(threshold, prevalence_frac, ~ {
    sum(!is.na(st$mean_abund) & st$mean_abund >= .x & st$prevalence >= .y &
        !is.na(st$cv) & st$cv >= MIN_CV)
  }), cv_filter = label_cv)

elbow_abund_combined <- bind_rows(elbow_abund_df, elbow_abund_cv_df)

p_elbow_abund <- ggplot(elbow_abund_combined,
                        aes(threshold, FeaturesRetained,
                            color = factor(prevalence_frac), linetype = cv_filter)) +
  geom_line(linewidth = 1.1) +
  scale_linetype_manual(values = setNames(c("solid", "dotted"), c(label_no_cv, label_cv)), name = NULL) +
  scale_x_log10() +
  geom_vline(xintercept = MIN_MEAN_ABUNDANCE, color = "black", linetype = "dashed") +
  theme_minimal() +
  labs(x = "Relative abundance threshold (log10 scale)",
       y = sprintf("Number of %s features passing filter", FEATURE_TABLE),
       color = "Min. prevalence\n(fraction of samples)",
       title = sprintf("Elbow plot for abundance filtering, with/without top %d%% CV filter",
                       MIN_CV_PERCENTILE))
ggsave(tag_filename("elbow_plot_abundance.pdf"), p_elbow_abund, width = 10, height = 6)
# ggsave(tag_filename("elbow_plot_abundance.png"), p_elbow_abund, width = 10, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Feature stats scatter (mean abundance vs. CV directly, both axes
# continuous) — a second, complementary view to the two elbow plots above:
# the elbow plots show retained *counts* as thresholds sweep, this shows
# every feature's actual position relative to both cutoffs at once.
# ---------------------------------------------------------------------------
write.csv(st, tag_filename("feature_stats.csv"), row.names = FALSE)

p_cv <- ggplot(st, aes(mean_abund, cv)) +
  geom_point(alpha = 0.35, size = 0.8, color = "steelblue") +
  scale_x_log10() +
  geom_vline(xintercept = MIN_MEAN_ABUNDANCE, color = "firebrick", linetype = "dashed") +
  geom_hline(yintercept = MIN_CV, color = "firebrick", linetype = "dashed") +
  theme_minimal() +
  labs(x = "Mean relative abundance (nonzero samples, log10 scale)",
       y = "Coefficient of variation (sd/mean, nonzero samples)",
       title = sprintf("%s: abundance vs. CV (n=%d features)", FEATURE_TABLE, nrow(st)))
ggsave(tag_filename("abundance_vs_cv.pdf"), p_cv, width = 8, height = 6)
# ggsave(tag_filename("abundance_vs_cv.png"), p_cv, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Abundance distribution histograms — permanent diagnostics for sanity-
# checking where MIN_MEAN_ABUNDANCE sits relative to the actual shape of the
# data (not just the elbow-plot summary of it). Two views:
#   1. Per-feature mean abundance (same values behind the elbow plot and the
#      abundance_vs_cv scatter) — shows whether the feature-level mean
#      abundances are unimodal (single unstable "steep region" around the
#      peak) or bimodal (a real, data-driven valley to threshold at).
#   2. All raw nonzero feature x sample values — can reveal multimodality
#      that per-feature averaging hides.
# ---------------------------------------------------------------------------
p_hist_mean <- ggplot(st, aes(mean_abund)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  scale_x_log10() +
  geom_vline(xintercept = MIN_MEAN_ABUNDANCE, color = "firebrick", linetype = "dashed") +
  theme_minimal() +
  labs(x = "Mean relative abundance (nonzero samples, log10 scale)",
       y = "Number of features",
       title = sprintf("%s: distribution of per-feature mean abundance", FEATURE_TABLE))
ggsave(tag_filename("abundance_histogram_mean.pdf"), p_hist_mean, width = 8, height = 6)
# ggsave(tag_filename("abundance_histogram_mean.png"), p_hist_mean, width = 8, height = 6, dpi = 300)

raw_vals <- mat[mat > 0]
p_hist_raw <- ggplot(data.frame(abund = raw_vals), aes(abund)) +
  geom_histogram(bins = 80, fill = "darkorange", color = "white") +
  scale_x_log10() +
  geom_vline(xintercept = MIN_MEAN_ABUNDANCE, color = "firebrick", linetype = "dashed") +
  theme_minimal() +
  labs(x = "Relative abundance (nonzero feature x sample values, log10 scale)",
       y = "Count",
       title = sprintf("%s: distribution of all nonzero abundance values", FEATURE_TABLE))
ggsave(tag_filename("abundance_histogram_raw.pdf"), p_hist_raw, width = 8, height = 6)
# ggsave(tag_filename("abundance_histogram_raw.png"), p_hist_raw, width = 8, height = 6, dpi = 300)
rm(raw_vals)  # can be large for genefamilies; drop once the plot is built

# ---------------------------------------------------------------------------
# Apply the filter
# ---------------------------------------------------------------------------
res <- filter_by_abundance_cv(mat, MIN_MEAN_ABUNDANCE, MIN_PREVALENCE, MIN_CV)
message(sprintf(
  "%d / %d features kept after abundance (>=%.2g), prevalence (>=%.2f), CV (>=%.2f) filters",
  nrow(res$mat), nrow(mat), MIN_MEAN_ABUNDANCE, MIN_PREVALENCE, MIN_CV))

if (nrow(res$mat) == 0) {
  stop("No features survived filtering — loosen MIN_MEAN_ABUNDANCE / MIN_CV / ",
       "MIN_PREVALENCE in config.R (see abundance_vs_cv.png).")
}

saveRDS(res$mat, FILTERED_FEATURES_RDS)
saveRDS(meta_aligned, META_ALIGNED_RDS)

writeLines(rownames(res$mat), tag_filename("surviving_pathway_ids_desc.txt"))
pathway_ids_clean <- trimws(sub(":.*$", "", rownames(res$mat)))
writeLines(pathway_ids_clean, tag_filename("surviving_pathway_ids_only.txt"))
message(sprintf("01_load_and_filter.R complete. Output in: %s", OUTPUT_DIR))