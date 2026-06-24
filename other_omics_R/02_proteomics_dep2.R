# =============================================================================
# 02_proteomics_dep2.R
#
# DEP2 differential-abundance workflow for proteomics / MS data. Unifies the
# original run_proteomics.R and run_proteomics_LC_GC.R via a loader switch:
#
#   PROTEOMICS$loader = "massspec"  per-sample long CSVs with Protein.Group,
#                                   Global.Q.Value, RawCount, contaminant;
#                                   runs the contaminant/Q-value QC.
#   PROTEOMICS$loader = "wide"      sample x feature CSVs (Sample column +
#                                   feature columns), merged on feature name;
#                                   QC is skipped.
#
# Downstream of the loader the two share one path: make_se -> filter_se ->
# normalize_vsn -> impute -> test_diff -> heatmap.
#
# Run:  Rscript 02_proteomics_dep2.R
# =============================================================================
source("config.R")
source(METADATA_UTILS)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(tidyr)
  library(ggplot2); library(tibble)
  library(DEP2); library(vsn)
  library(ComplexHeatmap); library(circlize)
})

cfg <- PROTEOMICS
res_dir <- cfg$output_dir
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Loaders: each returns a `combined_wide` data.frame with an ID column, a name
# column, and one column per sample (named like the metadata `label`), plus the
# integer `ecols` of the sample columns.
# ---------------------------------------------------------------------------
load_massspec <- function() {
  files <- list.files(file.path(cfg$parent_dir, "data"),
                      pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

  combined_df <- map_dfr(files, function(file) {
    df <- read_csv(file, col_types = cols(.default = "c"), quote = "\"") %>%
      select(-starts_with("..."))
    raw_col  <- grep("^[0-9]+_.*$|^Sample", colnames(df), value = TRUE)
    iBAQ_col <- grep("^iBAQ", colnames(df), value = TRUE)
    df <- df %>% rename(RawCount = all_of(raw_col),
                        iBAQ_NormCount = all_of(iBAQ_col))
    sample_name <- basename(file)
    for (s in cfg$sample_name_strip) sample_name <- gsub(s, "", sample_name, fixed = TRUE)
    df %>% mutate(SampleID = sample_name)
  })

  num_cols <- c("Global.Q.Value", "RawCount", "peptideCountsUnique",
                "peptideCountsAll", "iBAQ_NormCount", "nbTrypticPeptides")
  combined_df[num_cols] <- lapply(combined_df[num_cols], as.numeric)
  combined_df$Relative.abundance <- as.numeric(sub("%", "", combined_df$Relative.abundance))
  combined_df$contaminant <- as.logical(combined_df$contaminant)

  # --- QC: contaminant + Q-value bin barplot, then filter ---
  bins <- c(seq(0, 0.01, by = 0.001), Inf)
  qc <- combined_df %>%
    mutate(Category = ifelse(contaminant, "Contaminant", "Protein"),
           Qbin = ifelse(Category == "Protein",
                         as.character(cut(Global.Q.Value, breaks = bins,
                                          include.lowest = TRUE,
                                          labels = c(paste0(seq(0, 0.009, by = 0.001), "-",
                                                            seq(0.001, 0.01, by = 0.001)), ">0.01"))),
                         "Contaminant")) %>%
    count(Qbin)
  p <- ggplot(qc, aes(Qbin, n, fill = Qbin)) +
    geom_col(show.legend = FALSE) +
    labs(x = "Category / Q-value bin", y = "Number of rows",
         title = "Filtering Proteomics Data: Contaminants and Q-value bins") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "proteomics_qvalue_bins.pdf"), p, width = 8, height = 6)

  clean_df <- combined_df %>%
    filter(!contaminant, Global.Q.Value < cfg$qvalue_cutoff)

  wide <- clean_df %>%
    select(Protein.Group, Gene.names, SampleID, RawCount) %>%
    pivot_wider(names_from = SampleID, values_from = RawCount) %>%
    dplyr::rename(ID = Protein.Group, name = Gene.names)
  wide$ID <- paste(wide$ID, wide$name, sep = ";")
  list(wide = as.data.frame(wide), ecols = grep("SaMu", colnames(wide)))
}

load_wide <- function() {
  files <- list.files(file.path(cfg$parent_dir, "data"),
                      pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  wide <- map(files, function(file) {
    read_csv(file, show_col_types = FALSE) %>%
      column_to_rownames("Sample") %>% t() %>% as.data.frame() %>%
      rownames_to_column("name")
  }) %>%
    reduce(full_join, by = "name") %>%      # shared features merge; unique -> NA
    as.data.frame()
  wide$ID <- wide$name
  # Unify with metadata labels by removing underscores from sample column names.
  colnames(wide) <- gsub("_", "", colnames(wide))
  list(wide = wide, ecols = grep("^SaMu", colnames(wide)))
}

loaded <- switch(cfg$loader,
                 massspec = load_massspec(),
                 wide     = load_wide(),
                 stop("Unknown PROTEOMICS$loader: ", cfg$loader))
combined_wide <- loaded$wide
ecols <- loaded$ecols

# ---------------------------------------------------------------------------
# Metadata -> coldata aligned to the sample columns
# ---------------------------------------------------------------------------
meta_df <- load_samu_metadata(
  METADATA_CSV,
  keep_cols    = META_KEEP_COLS,
  fullsamu_col = FULLSAMU_COL,
  group_levels = GROUP_LEVELS,
  recode_smoke_alcohol = RECODE_SMOKE_ALCOHOL,
  add_label    = TRUE
)

common  <- intersect(colnames(combined_wide)[ecols], meta_df$label)
coldata <- meta_df[match(common, meta_df$label), ]
stopifnot(all(coldata$label %in% colnames(combined_wide)[ecols]))

# ---------------------------------------------------------------------------
# Build the SummarizedExperiment
# ---------------------------------------------------------------------------
coldata$condition <- as.character(factor(coldata[[GROUP_VAR]], levels = GROUP_LEVELS))
coldata$condition[is.na(coldata$condition)] <- "Unknown"
coldata$replicate <- as.numeric(gsub("SaMu", "", coldata$label))
coldata <- coldata[, c("label", "condition", "replicate")]

rownames(combined_wide) <- combined_wide$ID
combined_wide$name <- combined_wide$ID

se <- make_se(combined_wide, columns = ecols, expdesign = coldata)
se <- se[, colData(se)$condition != "Unknown"]

# ---------------------------------------------------------------------------
# Filtering threshold diagnostic + apply chosen filter
# ---------------------------------------------------------------------------
sweep_grid <- expand.grid(fraction = cfg$filter_sweep_fractions,
                          thr = cfg$filter_sweep_thresholds)
sweep_res <- sweep_grid %>%
  rowwise() %>%
  mutate(n_rows = nrow(filter_se(se, thr = thr, fraction = fraction))) %>%
  ungroup()
p <- ggplot(sweep_res, aes(fraction, n_rows, color = factor(thr))) +
  geom_line(size = 1) +
  labs(x = "missing fraction threshold (valid fraction required)",
       y = "features retained", color = "missing replicate threshold") +
  theme_minimal()
ggsave(file.path(res_dir, "dep2_se_filtering_thresholds.pdf"), p, width = 8, height = 6)

# NOTE: DEP2's `fraction` is the fraction of VALID values required (not missing).
se_filt <- filter_se(se, thr = cfg$filter_thr, fraction = cfg$filter_fraction)

pdf(file.path(res_dir, "dep2_filtering_effect.pdf"), width = 20, height = 6)
print(plot_frequency(se) + ggtitle("Identification overlap before filter"))
print(plot_frequency(se_filt) + ggtitle("Identification overlap after filter"))
dev.off()

# ---------------------------------------------------------------------------
# Normalize, inspect, impute
# ---------------------------------------------------------------------------
se_norm <- normalize_vsn(se_filt)
pdf(file.path(res_dir, "dep2_normalization_effect.pdf"), width = 6, height = 20)
print(plot_normalization(se_filt, se_norm)); dev.off()

pdf(file.path(res_dir, "meanSdPlot.pdf"), width = 8, height = 6)
meanSdPlot(assay(se_norm)); dev.off()

pdf(file.path(res_dir, "dep2_missing_values_heatmap.pdf"), width = 8, height = 6)
plot_missval(se_filt); dev.off()

se_imp <- DEP2::impute(se_norm, fun = cfg$impute_fun)

# Save imputed matrix (samples x features) and metadata.
proteomics_df <- as.data.frame(t(assay(se_imp)))
proteomics_df$Sample <- rownames(proteomics_df)
proteomics_df <- proteomics_df[, c("Sample", setdiff(colnames(proteomics_df), "Sample"))]
write.csv(proteomics_df,
          file.path(res_dir, paste0(cfg$experiment_name, "_dep2_vsn_imputed_matrix.csv")),
          row.names = FALSE)
write.csv(as.data.frame(colData(se_imp)),
          file.path(res_dir, paste0(cfg$experiment_name, "_dep2_vsn_imputed_metadata.csv")),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------
colData(se_imp)$condition <- relevel(factor(colData(se_imp)$condition),
                                     ref = GROUP_LEVELS[1])
diff <- test_diff(se_imp, type = "manual", test = cfg$de_test, fdr.type = "BH")
results <- add_rejections(diff)
write.csv(as.data.frame(rowData(results)),
          file.path(res_dir, paste0(cfg$experiment_name, "_dep2_results.csv")))

pdf(file.path(res_dir, paste0(cfg$experiment_name, "_dep2_volcano.pdf")),
    width = 8, height = 6)
print(plot_volcano(results, adjusted = FALSE, add_threshold_line = "intersect",
                   pCutoff = 0.05, fcCutoff = 1))
dev.off()

# ---------------------------------------------------------------------------
# Heatmap of top features by DE p-value
# ---------------------------------------------------------------------------
log_mat    <- log10(assay(se_filt) + 1e-6)
scaled_mat <- t(scale(t(log_mat)))

pcol <- paste0(cfg$de_test, "_p.val")
top_feats <- rowData(results) %>% as.data.frame() %>%
  arrange(.data[[pcol]]) %>% slice_head(n = cfg$heatmap_top_n) %>% pull(ID)
scaled_top <- scaled_mat[top_feats, ]

sample_anno <- data.frame(SarcStatus = colData(se_imp)$condition)
rownames(sample_anno) <- colnames(scaled_top)
ha <- HeatmapAnnotation(
  SarcStatus = sample_anno$SarcStatus,
  col = list(SarcStatus = c(setNames(c("steelblue", "orange"), GROUP_LEVELS),
                            Unknown = "gray")))

col_fun <- colorRamp2(c(min(scaled_top, na.rm = TRUE),
                        max(scaled_top, na.rm = TRUE)), c("black", "red"))

pdf(file.path(res_dir, paste0(cfg$experiment_name, "_heatmap.pdf")),
    width = 15, height = 12)
draw(Heatmap(scaled_top, col = col_fun, na_col = "gray80",
             show_row_names = TRUE, row_names_side = "left",
             row_names_gp = gpar(fontsize = 10),
             name = "scaled Log10(Intensity)", top_annotation = ha,
             show_column_names = TRUE, cluster_columns = FALSE,
             column_order = order(colnames(scaled_top)),
             clustering_distance_rows = "euclidean",
             clustering_method_rows = "complete",
             column_title = "Samples", row_title = "Features"))
dev.off()

message("02_proteomics_dep2.R complete (loader: ", cfg$loader, ").")
