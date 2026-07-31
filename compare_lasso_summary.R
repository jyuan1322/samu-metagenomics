# =============================================================================
# compare_lasso_summary.R
#
# Compare a chosen parameter (e.g. odds_ratio, nonzero_fraction, coefficient)
# between two lasso summary tables (e.g. two different contrasts, outcome
# encodings, or covariate adjustments run on the same features) and plot a
# scatter of param_1 vs param_2, joined by feature ID.
#
# Unlike compare_deseq2_lfc.R, the compared column is a CONFIG parameter
# rather than hardcoded — set compare_col to whatever column you want to
# compare (odds_ratio, nonzero_fraction, coefficient, ...).
#
# ID column handling: if id_col (below) names an existing column, that's
# used directly. Otherwise, falls back to auto-detecting an unnamed first
# (row-name-style) column — the shape you get from write.csv(..., row.names
# = TRUE), which several of your other scripts already produce output in.
#
# Edit CONFIG below, then:  Rscript compare_lasso_summary.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2)
})

# ---------------------------------------------------------------------------
CONFIG <- list(
  file1  = "/data/local/jy1008/SaMu/results/latest/metagenomics_ml_test2/lasso_L1_stratifiedkfold_i111_o222_coef_summary.csv",
  file2  = "/data/local/jy1008/SaMu/results/latest/metagenomics_ml_extreme_cases_test2/lasso_L1_stratifiedkfold_i111_o222_coef_summary.csv",
  label1 = "Model 1: full data",
  label2 = "Model 2: severe cases only",

  output_dir    = "/data/local/jy1008/SaMu/results/latest",
  output_prefix = "lasso_compare",

  # Delimiter for both input files: "auto" sniffs from the first line
  # (counts tabs vs commas), or set explicitly ("\t" or ",").
  delim = "auto",

  # Name of the feature/taxon identifier column. If this doesn't exist in
  # the file, falls back to auto-detecting an unnamed first column.
  id_col = "feature",

  # The column to compare/plot. Set to whatever your summary table calls
  # it: "odds_ratio", "nonzero_fraction", "coefficient", etc.
  compare_col = "odds_ratio",

  # Optional highlighting column + threshold, analogous to padj for DESeq2
  # (e.g. a selection/nonzero frequency across CV folds). Set highlight_col
  # to NULL to skip highlighting and plot a single-color scatter instead.
  # highlight_col       = "nonzero_fraction",
  highlight_col       = NULL,
  highlight_threshold = 0.5,
  # "above": highlight_col > threshold counts as "selected".
  # "below": highlight_col < threshold counts as "selected" (e.g. for a
  #          p-value-like column where smaller = more confident).
  highlight_direction = "above"
)
# ---------------------------------------------------------------------------

dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# sniff_delim — guess comma vs tab from the first line of the file.
# ---------------------------------------------------------------------------
sniff_delim <- function(path) {
  first_line <- readLines(path, n = 1)
  n_tab   <- lengths(regmatches(first_line, gregexpr("\t", first_line)))
  n_comma <- lengths(regmatches(first_line, gregexpr(",", first_line)))
  if (n_tab > n_comma) "\t" else ","
}

# ---------------------------------------------------------------------------
# read_summary_csv — read a summary table, resolve the id column to
# "feature" regardless of whether it's an explicitly named column or an
# R-auto-detected unnamed first (row-name-style) column.
# ---------------------------------------------------------------------------
read_summary_csv <- function(path, delim = "auto", id_col = "feature") {
  if (delim == "auto") delim <- sniff_delim(path)
  df <- read.delim(path, sep = delim, header = TRUE, check.names = FALSE)

  if (id_col %in% colnames(df)) {
    df <- df %>% rename(feature = all_of(id_col))
  } else {
    auto_rownames <- !all(rownames(df) == as.character(seq_len(nrow(df))))
    if (auto_rownames) {
      df <- df %>% rownames_to_column("feature")
    } else {
      stop(sprintf(
        "id_col '%s' not found in %s, and no unnamed row-name-style column detected.\nColumns found: %s",
        id_col, path, paste(colnames(df), collapse = ", ")))
    }
  }
  df
}

d1 <- read_summary_csv(CONFIG$file1, CONFIG$delim, CONFIG$id_col)
d2 <- read_summary_csv(CONFIG$file2, CONFIG$delim, CONFIG$id_col)

for (df_name in c("d1", "d2")) {
  df <- get(df_name)
  if (!CONFIG$compare_col %in% colnames(df))
    stop(sprintf("compare_col '%s' not found in %s. Columns found: %s",
                 CONFIG$compare_col, df_name, paste(colnames(df), collapse = ", ")))
}

message(sprintf("File 1 (%s): %d features", CONFIG$label1, nrow(d1)))
message(sprintf("File 2 (%s): %d features", CONFIG$label2, nrow(d2)))

# ---------------------------------------------------------------------------
# Inner join on shared features only — comparison isn't meaningful for
# features present in just one file (counts are reported, not silently
# discarded).
# ---------------------------------------------------------------------------
merged <- inner_join(d1, d2, by = "feature", suffix = c("_1", "_2"))

n_only1 <- sum(!d1$feature %in% d2$feature)
n_only2 <- sum(!d2$feature %in% d1$feature)
message(sprintf("Shared features: %d | only in file 1: %d | only in file 2: %d",
                nrow(merged), n_only1, n_only2))

write.csv(merged,
          file.path(CONFIG$output_dir, paste0(CONFIG$output_prefix, "_merged.csv")),
          row.names = FALSE)

compare_col_1 <- paste0(CONFIG$compare_col, "_1")
compare_col_2 <- paste0(CONFIG$compare_col, "_2")

# ---------------------------------------------------------------------------
# Optional highlighting flag + correlation for the plot annotation.
# ---------------------------------------------------------------------------
if (!is.null(CONFIG$highlight_col)) {
  hl_col_1 <- paste0(CONFIG$highlight_col, "_1")
  hl_col_2 <- paste0(CONFIG$highlight_col, "_2")

  is_selected <- function(x) {
    if (CONFIG$highlight_direction == "above") x > CONFIG$highlight_threshold
    else x < CONFIG$highlight_threshold
  }

  hl_both  <- "Selected in both"
  hl_1only <- paste0("Selected in ", CONFIG$label1, " only")
  hl_2only <- paste0("Selected in ", CONFIG$label2, " only")
  hl_none  <- "Not selected"

  merged <- merged %>%
    mutate(sig = case_when(
      is_selected(.data[[hl_col_1]]) & is_selected(.data[[hl_col_2]]) ~ hl_both,
      is_selected(.data[[hl_col_1]]) ~ hl_1only,
      is_selected(.data[[hl_col_2]]) ~ hl_2only,
      TRUE ~ hl_none
    ))

  sig_colors <- setNames(
    c("firebrick", "dodgerblue3", "darkorange", "grey70"),
    c(hl_both, hl_1only, hl_2only, hl_none)
  )
} else {
  merged$sig <- "all"
  sig_colors <- c("all" = "grey30")
}

cor_val <- cor(merged[[compare_col_1]], merged[[compare_col_2]],
              method = "spearman", use = "complete.obs")

# ---------------------------------------------------------------------------
# Scatterplot: compare_col (file1) vs compare_col (file2)
# ---------------------------------------------------------------------------
lim <- max(abs(c(merged[[compare_col_1]], merged[[compare_col_2]])), na.rm = TRUE) * 1.05

p <- ggplot(merged, aes(.data[[compare_col_1]], .data[[compare_col_2]], color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = sig_colors) +
  coord_equal(xlim = c(0.75, 1.1), ylim = c(0.75, 1.1)) +
  theme_classic() +
  labs(
    title = sprintf("%s comparison: %s vs %s", CONFIG$compare_col, CONFIG$label1, CONFIG$label2),
    subtitle = sprintf("Spearman correlation = %.3f (n = %d shared features)",
                       cor_val, nrow(merged)),
    x = paste0(CONFIG$label1, " ", CONFIG$compare_col),
    y = paste0(CONFIG$label2, " ", CONFIG$compare_col),
    color = NULL
  ) +
  { if (is.null(CONFIG$highlight_col)) guides(color = "none") }

out_base <- file.path(CONFIG$output_dir, paste0(CONFIG$output_prefix, "_", CONFIG$compare_col, "_scatter"))
ggsave(paste0(out_base, ".pdf"), p, width = 7, height = 6.5)
ggsave(paste0(out_base, ".png"), p, width = 7, height = 6.5, dpi = 300)

message("compare_lasso_summary.R complete.")