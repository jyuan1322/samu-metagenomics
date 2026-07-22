# =============================================================================
# 01b_quorum.R
#
# Quorum-sensing-molecule analysis: join quorum assay table to metadata,
# group by QSP / Species / Microbial_target, log-transform and z-score
# scale raw (summed) quant values, and produce a clustered heatmap per
# grouping level. Config in config.R under NMR. Runs independently of NMR
# data — does not require nmr_csv to be present.
#
# Run:  Rscript 01b_quorum.R
# =============================================================================
source("config.R")
source(METADATA_UTILS)
source("matrix_utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr); library(ggplot2)
  library(pheatmap); library(tibble)
})

cfg <- NMR
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(cfg$output_dir)

# ---------------------------------------------------------------------------
# Load metadata + quorum assay table, build a merged frame
# ---------------------------------------------------------------------------
# Metadata is recoded from the raw CSV (not the metagenomics output) so that
# quorum samples outside the metagenomics cohort are still covered.
meta_df <- load_samu_metadata(
  METADATA_CSV,
  keep_cols    = META_KEEP_COLS,
  fullsamu_col = FULLSAMU_COL,
  group_levels = GROUP_LEVELS,
  recode_smoke_alcohol = RECODE_SMOKE_ALCOHOL
)

quorum_concs <- read.csv(file.path(cfg$input_dir, cfg$quorum_csv))
quorum_cols  <- str_subset(colnames(quorum_concs), cfg$quorum_col_regex)

# Below-detection-limit placeholder -> 0 in quorum columns.
quorum_concs[quorum_cols] <- lapply(quorum_concs[quorum_cols], function(x) {
  x[x == cfg$quorum_below_detection] <- "0"; x
})

# load_samu_metadata already applied the Full.SaMu filter; a left join keeps
# every retained metadata sample and attaches quorum values where present.
merged_df <- meta_df %>%
  merge(quorum_concs[, c("record_id", quorum_cols)], by = "record_id", all.x = TRUE)
merged_df[[GROUP_VAR]] <- factor(merged_df[[GROUP_VAR]], levels = GROUP_LEVELS)

write.csv(merged_df, "samu_quorum_FullSaMu.csv", row.names = FALSE, quote = TRUE)

# ---------------------------------------------------------------------------
# QSP -> Species / Microbial_target grouping reference.
#
# Character columns are trimmed (the source table has stray trailing spaces,
# e.g. "Streptococcus mitis " vs "Streptococcus mitis", which would otherwise
# split one species into two groups). Blank Species/Microbial_target entries
# become NA so they can be dropped from the relevant grouping level below
# without affecting the other levels or the raw QSP-level analysis.
# ---------------------------------------------------------------------------
grouping_df <- read.csv(file.path(cfg$input_dir, cfg$quorum_grouping_file),
                         stringsAsFactors = FALSE, check.names = FALSE)
colnames(grouping_df) <- trimws(colnames(grouping_df))

required_cols <- c("QSP", "Species", "Microbial_target")
missing_cols  <- setdiff(required_cols, colnames(grouping_df))
if (length(missing_cols) > 0) {
  stop("quorum_grouping_file is missing expected column(s): ",
       paste(missing_cols, collapse = ", "),
       "\nColumns found instead: ", paste(colnames(grouping_df), collapse = ", "))
}

grouping_df[] <- lapply(grouping_df, function(x) if (is.character(x)) trimws(x) else x)
grouping_df$Species[grouping_df$Species == ""] <- NA
grouping_df$Microbial_target[grouping_df$Microbial_target == ""] <- NA

# ---------------------------------------------------------------------------
# aggregate_quorum_cols — collapse QSP quant columns to a grouping level by
# summing (quant*prob) across QSP columns that share a group value. QSP
# columns with no assignment at this level (blank in grouping_df, or absent
# from it entirely) are dropped from this level only.
#
# Args:
#   mat         : sample x QSP matrix, colnames like "<QSP>_Quant.Prob"
#   col_regex   : regex used to strip the QSP suffix off colnames
#   grouping_df : QSP -> Species / Microbial_target reference table
#   group_col   : "QSP" (pass-through, no grouping file needed),
#                 "Species", or "Microbial_target"
# ---------------------------------------------------------------------------
aggregate_quorum_cols <- function(mat, col_regex, grouping_df, group_col) {
  if (group_col == "QSP") return(mat)

  qsp_id    <- sub(col_regex, "", colnames(mat))
  group_val <- grouping_df[[group_col]][match(qsp_id, grouping_df$QSP)]

  keep <- !is.na(group_val)
  if (!all(keep)) {
    message(sprintf("  dropping %d/%d QSP columns with no %s assignment: %s",
                     sum(!keep), length(keep), group_col,
                     paste(qsp_id[!keep], collapse = ", ")))
  }
  mat       <- mat[, keep, drop = FALSE]
  group_val <- group_val[keep]

  agg <- sapply(split(seq_len(ncol(mat)), group_val), function(idx) {
    rowSums(mat[, idx, drop = FALSE])
  })
  rownames(agg) <- rownames(mat)
  agg
}

# ---------------------------------------------------------------------------
# Quorum sensing: log-transform raw (grouped) values directly, matching the
# NMR approach — no compositional/relative-abundance normalization.
# Run once per grouping level (QSP columns as-is, summed by Species, summed
# by Microbial_target); each level gets its own CSV and heatmap.
# ---------------------------------------------------------------------------
quorum_mat <- build_matrix(merged_df, quorum_cols)
quorum_mat[is.na(quorum_mat)] <- 0                       # unique to quorum data
quorum_mat <- quorum_mat[rowSums(quorum_mat) > 0, , drop = FALSE]
quorum_mat <- quorum_mat[, colSums(quorum_mat) > 0, drop = FALSE]

heat_colors <- colorRampPalette(c("black", "red"))(100)

grouping_levels <- c(qsp = "QSP", species = "Species",
                     microbial_target = "Microbial_target")

for (level_slug in names(grouping_levels)) {
  group_col <- grouping_levels[[level_slug]]
  message("Quorum grouping level: ", group_col)

  level_mat <- aggregate_quorum_cols(quorum_mat, cfg$quorum_col_regex,
                                     grouping_df, group_col)
  level_mat <- level_mat[, colSums(level_mat) > 0, drop = FALSE]
  level_mat <- level_mat[rowSums(level_mat) > 0, , drop = FALSE]  # a sample
  # can end up all-zero at this level even if not at QSP level, e.g. if its
  # only signal was on QSPs excluded from this grouping.

  # NOTE: quorum values (quant*prob, summed within a group where applicable)
  # are log-transformed directly, the same way NMR concentrations are —
  # NOT converted to relative abundances first. No compositional constraint
  # is assumed for this assay.
  log_quorum <- log10(level_mat + cfg$epsilon)
  log_quorum_scaled_full <- scale(log_quorum, center = TRUE, scale = TRUE)

  gaq <- group_annotation(log_quorum_scaled_full, merged_df, GROUP_VAR, GROUP_LEVELS,
                          cluster_within = TRUE)
  log_quorum_scaled <- log_quorum_scaled_full[gaq$order, , drop = FALSE]

  write.csv(log_quorum_scaled, paste0("samu_quorum_log_scaled_", level_slug, ".csv"),
            row.names = TRUE, quote = TRUE)

  # hard-coded colors
  ann_colors <- list(sarc_status_bin = setNames(c("#ff9289ff", "#00dae0ff"), GROUP_LEVELS))
  pheatmap(log_quorum_scaled, color = heat_colors, scale = "none",
           clustering_method = "ward.D2", cluster_cols = TRUE, cluster_rows = FALSE,
           annotation_row = gaq$anno, gaps_row = gaq$gap, fontsize = 10,
           border_color = NA, angle_col = 45,
           annotation_colors = ann_colors,
           filename = paste0("samu_quorum_log_scaled_", level_slug, ".pdf"),
           width = 12, height = 12)

  # -------------------------------------------------------------------------
  # PCA (scaled log values), colored by group / age / total signal.
  # total_signal is the quorum analogue of NMR's "read depth": rowSums of
  # the raw (pre-log) grouped values, i.e. total quant*prob per sample.
  # -------------------------------------------------------------------------
  log_scaled_clean <- log_quorum_scaled[, colSums(is.na(log_quorum_scaled)) == 0,
                                        drop = FALSE]
  pca_res <- prcomp(log_scaled_clean, center = FALSE, scale. = FALSE)
  pca_df  <- as.data.frame(pca_res$x)
  pca_df$Sample <- rownames(pca_df)
  pca_df <- cbind(pca_df,
                  merged_df[match(pca_df$Sample, merged_df$record_id),
                            c(GROUP_VAR, "age_def"), drop = FALSE])
  pca_df$age_def <- as.numeric(pca_df$age_def)  # metadata CSV can load this as character
  pca_df$total_signal <- rowSums(level_mat)[pca_df$Sample]

  ggsave(paste0("samu_quorum_pca_1_color_by_sarc_status_", level_slug, ".pdf"),
         ggplot(pca_df, aes(PC1, PC2, color = .data[[GROUP_VAR]])) +
           geom_point(size = 3) + theme_classic() +
           labs(title = paste0("PCA colored by group (", group_col, " level)"),
                color = "Group"),
         width = 6, height = 5)
  ggsave(paste0("samu_quorum_pca_color_by_age_", level_slug, ".pdf"),
         ggplot(pca_df, aes(PC1, PC2, color = age_def)) +
           geom_point(size = 3) + scale_color_viridis_c() + theme_classic() +
           labs(title = paste0("PCA colored by age (", group_col, " level)"),
                color = "Age (years)"),
         width = 6, height = 5)
  ggsave(paste0("samu_quorum_pca_color_by_total_signal_", level_slug, ".pdf"),
         ggplot(pca_df, aes(PC1, PC2, color = total_signal)) +
           geom_point(size = 3) + scale_color_viridis_c(option = "C") +
           theme_classic() +
           labs(title = paste0("PCA colored by total quant*prob signal (", group_col, " level)"),
                color = "Total signal"),
         width = 6, height = 5)

  # -------------------------------------------------------------------------
  # Per-feature Wilcoxon tests + boxplots (raw, un-scaled grouped values).
  # "feature" is a QSP id, species, or microbial target depending on level.
  # -------------------------------------------------------------------------
  long_df <- as.data.frame(level_mat) %>%
    rownames_to_column("record_id") %>%
    mutate(record_id = as.character(record_id)) %>%
    left_join(merged_df %>%
                mutate(record_id = as.character(record_id)) %>%
                select(record_id, all_of(GROUP_VAR)),
              by = "record_id") %>%
    pivot_longer(-c(record_id, all_of(GROUP_VAR)),
                 names_to = "feature", values_to = "value")

  pvals <- long_df %>%
    group_by(feature) %>%
    summarise(p_value = wilcox.test(value ~ .data[[GROUP_VAR]])$p.value,
              .groups = "drop") %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))

  y_positions <- long_df %>%
    mutate(log_value = log10(value + 0.01)) %>%
    group_by(feature) %>%
    summarise(y_pos = max(log_value, na.rm = TRUE) * 1.1, .groups = "drop")

  pvals_plot <- pvals %>%
    left_join(y_positions, by = "feature") %>%
    mutate(label = paste0("p = ", signif(p_value, 2),
                          "\n(padj = ", signif(p_adj, 2), ")"))

  # Order features by raw p-value on the x-axis.
  feat_order <- pvals_plot %>% arrange(p_value) %>% pull(feature)
  long_df$feature    <- factor(long_df$feature, levels = feat_order)
  pvals_plot$feature <- factor(pvals_plot$feature, levels = feat_order)

  n_feat <- n_distinct(long_df$feature)
  p <- ggplot(long_df, aes(feature, log10(value + cfg$epsilon), fill = .data[[GROUP_VAR]])) +
    geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
    geom_jitter(aes(color = .data[[GROUP_VAR]]),
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                alpha = 0.5, size = 0.7, show.legend = FALSE) +
    geom_text(data = pvals_plot, aes(feature, y_pos, label = label),
              inherit.aes = FALSE, size = 3.5) +
    theme_classic() +
    labs(x = group_col, y = paste0("log10(Quant*Prob + ", cfg$epsilon, ")"), fill = "Group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("quorum_wilcoxon_boxplot_", level_slug, ".pdf"), p,
         width = max(6, n_feat * 1.5), height = 6, dpi = 300)
}

message("01b_quorum.R complete.")