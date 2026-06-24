library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)
library(tibble)

meta_df <- read.csv("/data/local/jy1008/SaMu/results/latest/metagenomics_R/SaMu_sarcopeniestatus_majorcovariates_v13_16012026_FullSaMu_only.csv",
                header = TRUE, stringsAsFactors = FALSE)
meta_df$label <- paste0("SaMu", meta_df$record_id)

# Filter metadata for important covariates
meta_df <- meta_df[, c(
  "record_id",
  "label",
  "Block",
  "Full.SaMu",
  "sarc_status",
  "EWGSOP_strength",
  "EWGSOP_mass",
  "EWGSOP_performance",
  "age_def",
  "sex",
  "smke",
  "alco",
  "nutr_score",
  "bmi",
  "stvol"
)]
meta_df <- meta_df %>% filter(Full.SaMu == 1)
# NOTE: 5/8/2026 - smoke and alcohol variables are already binarized, so no need to convert.
# meta_df$smke <- ifelse(meta_df$smke %in% c("1", "2", "3"), 1,
#                        ifelse(meta_df$smke == "0", 0, NA))
# meta_df$alco <- ifelse(meta_df$alco %in% c("1", "2"), 0,
#                        ifelse(meta_df$alco %in% c("3", "4"), 1, NA))
print("Converting sarc_status to numeric, converting non-numeric values to NA")
meta_df$sarc_status <- as.numeric(meta_df$sarc_status)
# binarize sarc_status
meta_df$sarc_status_bin <- ifelse(meta_df$sarc_status > 0, "Sarc", "NoSarc")
meta_df$sarc_status_bin <- factor(meta_df$sarc_status_bin, levels = c("NoSarc", "Sarc"))




# parent_dir <- "/data/local/jy1008/SaMu/proteomics/GC_MS"
# output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_GC-MS"
# experiment_name <- "GC_MS"
# parent_dir <- "/data/local/jy1008/SaMu/proteomics/LC_MS_neg"
# output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_neg/clustered_all"
# experiment_name <- "LC_MS_neg_clustered"
# output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_neg/clustered_1567"
# experiment_name <- "LC_MS_neg_1567"
# output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_neg/clustered_234"
# experiment_name <- "LC_MS_neg_234"
# output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_neg/clustered_89"
# experiment_name <- "LC_MS_neg_89"

parent_dir <- "/data/local/jy1008/SaMu/proteomics/LC_MS_pos"
# output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_pos/clustered_all"
# experiment_name <- "LC_MS_pos_clustered"
# output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_pos/clustered_1567"
# experiment_name <- "LC_MS_pos_1567"
# output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_pos/clustered_234"
# experiment_name <- "LC_MS_pos_234"
output_dir <- "/data/local/jy1008/SaMu/results/latest/proteomics_LC-MS_pos/clustered_89"
experiment_name <- "LC_MS_pos_89"


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# files <- list.files(file.path(parent_dir, "data"),
#                     pattern = "\\.csv$",
#                     recursive = TRUE,
#                     full.names = TRUE)
# print(files)

# for LC-MS test, single batch
# files <- c("/data/local/jy1008/SaMu/proteomics/LC_MS_neg/data/LCMS_processed_blocks1,5,6,7_SaMu_2_28and165_297_neg.csv")
# files <- c("/data/local/jy1008/SaMu/proteomics/LC_MS_neg/data/LCMS_processed_blocks2,3,4_SaMu_29_164_neg.csv")
# files <- c("/data/local/jy1008/SaMu/proteomics/LC_MS_neg/data/LCMS_processed_blocks8,9,1PY2_SaMu_298_371and5PY2_25PY2_neg.csv")

# files <- c("/data/local/jy1008/SaMu/proteomics/LC_MS_pos/data/LCMS_processed_blocks1,5,6,7_SaMu_2_28and165_297_pos.csv")
# files <- c("/data/local/jy1008/SaMu/proteomics/LC_MS_pos/data/LCMS_processed_blocks2,3,4_SaMu_29_164_pos.csv")
files <- c("/data/local/jy1008/SaMu/proteomics/LC_MS_pos/data/LCMS_processed_blocks8,9,1PY2_SaMu_298_371and5PY2_25PY2_pos.csv")

# ----------------------------------------
# Table: which features are in which files
# ----------------------------------------
feature_file_map <- map_dfr(files, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  
  # The non-Sample columns are the features/metabolites
  features <- setdiff(colnames(df), "Sample")
  
  tibble(
    Metabolite = features,
    File = basename(file)
  )
})

# Pivot to wide: rows = features, cols = files, values = presence (0/1)
feature_presence_table <- feature_file_map %>%
  mutate(Present = 1) %>%
  pivot_wider(
    names_from = File,
    values_from = Present,
    values_fill = 0
  )

# Summary: how many features per file
file_summary <- feature_file_map %>%
  count(File, name = "n_features")
print(file_summary)

# Summary: how many files each feature appears in
feature_summary <- feature_presence_table %>%
  mutate(n_files = rowSums(select(., -Metabolite))) %>%
  count(n_files, name = "n_features")
print(feature_summary)

# Save the table
write_csv(feature_presence_table, file.path(output_dir, paste0(experiment_name, "_feature_file_presence.csv")))




# ----------------------------------------
# Combine files: shared features get one row, batch-unique get NAs
# ----------------------------------------
combined_wide <- map(files, function(file) {
  df <- as.data.frame(read_csv(file, show_col_types = FALSE))
  rownames(df) <- df$Sample
  df$Sample <- NULL
  df_t <- as.data.frame(t(df))
  df_t$name <- rownames(df_t)
  df_t
}) %>%
  reduce(full_join, by = "name")

combined_wide <- as.data.frame(combined_wide)

# unify with meta_df labels by removing underscores from column names
colnames(combined_wide) <- gsub("_", "", colnames(combined_wide))
ecols <- grep("^SaMu", colnames(combined_wide))


# combined_wide <- combined_df %>%
#   select(Protein.Group, SampleID, RawCount) %>%
#   tidyr::pivot_wider(names_from = SampleID,
#                      values_from = RawCount)

# ----------------------------------------
# QC and visualization
# ----------------------------------------

# # Make numeric Q-values
# combined_df <- combined_df %>%
#   mutate(
#     Qval = as.numeric(Global.Q.Value),
#     Category = ifelse(contaminant, "Contaminant", "Protein")
#   )

# # Define Q-value bins for non-contaminants
# # bins up to 0.01, then one bin for everything above
# bins <- c(seq(0, 0.01, by = 0.001), Inf)  # last bin = >0.01

# combined_df <- combined_df %>%
#   mutate(
#     Qbin = case_when(
#       Category == "Protein" ~ cut(
#         Qval,
#         breaks = bins,
#         include.lowest = TRUE,
#         labels = c(paste0(seq(0, 0.009, by = 0.001), "-", seq(0.001, 0.01, by = 0.001)), ">0.01")
#       ),
#       TRUE ~ "Contaminant"
#     )
#   )

# # Count rows per category/bin
# plot_data <- combined_df %>%
#   group_by(Qbin) %>%
#   summarise(count = n()) %>%
#   ungroup()

# p <- ggplot(plot_data, aes(x = Qbin, y = count, fill = Qbin)) +
#   geom_col(show.legend = FALSE) +
#   scale_fill_manual(values = c("Contaminant" = "red", rep("steelblue", length(unique(plot_data$Qbin)) - 1))) +
#   labs(
#     x = "Category / Q-value bin",
#     y = "Number of rows",
#     title = "Filtering Proteomics Data: Contaminants and Q-value bins"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(file.path(parent_dir, "results", "proteomics_qvalue_bins.pdf"), p, width = 8, height = 6)


# Perform the actual QC

# # filter out contaminants and keep only those with Global.Q.Value < 0.01
# clean_df <- combined_df %>%
#   filter(!contaminant) %>% 
#   filter(as.numeric(Global.Q.Value) < 0.01)

# combined_wide <- clean_df %>%
#   select(Protein.Group, Gene.names, SampleID, RawCount) %>%  # use RawCount
#   pivot_wider(
#     names_from = SampleID,
#     values_from = RawCount
#   )

# # combined_wide should have Protein.Group and Protein.Names
# combined_wide <- combined_wide %>%
#   dplyr::rename(ID = Protein.Group,
#                 name = Gene.names)  # or Gene.names if you prefer

# now make unique
# combined_unique <- make_unique(
#   combined_wide,
#   "ID",
#   "name"
# )
# combined_wide$ID <- paste(combined_wide$ID, combined_wide$name, sep = ";")
# 
# ecols <- grep("SaMu", colnames(combined_wide))

# ----------------------------------------
# Metadata import
# ----------------------------------------


# intensity_matrix: columns are sample IDs
# coldata: metadata table

# Keep only samples that exist in intensity_matrix
coldata <- meta_df %>%
  filter(label %in% colnames(combined_wide)[ecols])

common <- intersect(colnames(combined_wide)[ecols], coldata$label)

coldata <- coldata[match(common, coldata$label), ]



# # Optional: reorder rows to match the column order in intensity_matrix
# coldata <- coldata[match(colnames(combined_wide)[ecols], coldata$label), ]

# # check that sample names match intensity matrix
all(coldata$label %in% colnames(combined_wide)[ecols])

# import metadata
# coldata <- data.frame(
#   label = colnames(intensity_matrix),
#   condition = c("Control","Control","Treatment","Treatment"), # replace with your design
#   replicate = c(1,2,1,2)
# )




# calculate the correlation of values
rownames(combined_wide) <- combined_wide$name
combined_wide$name <- NULL
# replace NA with 0
combined_wide[is.na(combined_wide)] <- 0


corr_thresh = 0.9
mat <- as.matrix(combined_wide)

# only keep peaks that are present in at least 25% of samples (to reduce noise from very sparse peaks)
peak_filter_threshold <- 0.25  # 25%
keep <- rowSums(mat != 0) / ncol(mat) >= peak_filter_threshold
mat <- mat[keep, ]

cor_mat <- cor(t(mat),
               method = "pearson",
               use = "pairwise.complete.obs")
rownames(cor_mat) <- rownames(mat)
colnames(cor_mat) <- rownames(mat)

library(pheatmap)

# Zero out the diagonal for display purposes
cor_mat_display <- cor_mat
diag(cor_mat_display) <- 0

pdf(file.path(output_dir, "lc_ms_neg_correlation_heatmap.pdf"), width = 20, height = 20)
pheatmap(cor_mat_display,
         clustering_distance_rows = as.dist(1 - cor_mat),
         clustering_distance_cols = as.dist(1 - cor_mat),
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()
pdf(file.path(output_dir, "lc_ms_neg_correlation_heatmap_highlight.pdf"), width = 20, height = 20)
cor_highlight <- cor_mat
cor_highlight[cor_highlight <= corr_thresh] <- NA

dist_mat <- as.dist(1 - cor_mat)

diag(cor_highlight) <- NA
pheatmap(cor_highlight,
         clustering_distance_rows = dist_mat,
         clustering_distance_cols = dist_mat,
         na_col = "white",
         color = colorRampPalette(c("yellow", "red"))(100),
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()


adj_w <- cor_mat
adj_w[adj_w < corr_thresh] <- 0
diag(adj_w) <- 0

# build graph and cluster peaks
library(igraph)

g <- graph_from_adjacency_matrix(adj_w,
                                 mode = "undirected",
                                 weighted = TRUE,
                                 diag = FALSE)

clusters <- components(g)

peak_clusters <- data.frame(
  peak = names(clusters$membership),
  cluster = clusters$membership
)

# Keep only clusters with >1 member
peak_clusters_filtered <- peak_clusters[
  clusters$csize[peak_clusters$cluster] > 1,
]

included_peaks <- peak_clusters_filtered$peak

# ----------------------------------------
# Summary: samples, features, clusters
# ----------------------------------------
n_samples       <- ncol(mat)
n_features_raw  <- nrow(mat)
n_clusters      <- length(unique(peak_clusters_filtered$cluster))
n_singleton     <- sum(clusters$csize == 1)  # peaks not in any cluster

sink(file.path(output_dir, paste0(experiment_name, "_pipeline_summary.txt")))
cat("=== Pipeline Summary ===\n")
cat(sprintf("  Samples:                        %d\n", n_samples))
cat(sprintf("  Raw features (post-filter):     %d\n", n_features_raw))
cat(sprintf("  Features in clusters (>1 peak): %d\n", length(included_peaks)))
cat(sprintf("  Singleton peaks (excluded):     %d\n", n_singleton))
cat(sprintf("  Clusters generated:             %d\n", n_clusters))
cat("========================\n")
sink()

cor_subset <- cor_mat[included_peaks, included_peaks]
pdf(file.path(output_dir, "lc_ms_neg_correlation_heatmap_highlight_filtered.pdf"), width = 20, height = 20)
cor_highlight <- cor_subset
cor_highlight[cor_highlight <= corr_thresh] <- NA

dist_mat <- as.dist(1 - cor_subset)

diag(cor_highlight) <- NA
pheatmap(cor_highlight,
         clustering_distance_rows = dist_mat,
         clustering_distance_cols = dist_mat,
         na_col = "white",
         color = colorRampPalette(c("yellow", "red"))(100),
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

# given peak clusters, reformat combined_wide to have one row per cluster,
# taking the mean intensity in each group as representative of the cluster
library(dplyr)
library(tidyr)

df_long <- as.data.frame(mat) %>%
  mutate(peak = rownames(mat)) %>%
  pivot_longer(
    cols = -peak,
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  inner_join(peak_clusters_filtered, by = "peak")

# Max intensity per cluster/sample
# df_clustered <- df_long %>%
#   group_by(cluster, Sample) %>%
#   summarise(Intensity = max(Intensity, na.rm = TRUE), .groups = "drop")
df_clustered <- df_long %>%
  group_by(cluster, Sample) %>%
  summarise(Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")


df_wide <- df_clustered %>%
  pivot_wider(
    names_from = Sample,
    values_from = Intensity
  )
colnames(df_wide)[1] <- "name"
df_wide["name"] <- paste0("Cluster", df_wide$name)
combined_wide <- df_wide


# add the "ID" column for DEP2 compatibility
combined_wide <- combined_wide %>%
  mutate(ID = name) %>%
  select(ID, everything()) # moves ID to the front
ecols <- grep("^SaMu", colnames(combined_wide))

# ----------------------------------------
# Metadata import (again, after clustering)
# ----------------------------------------


# intensity_matrix: columns are sample IDs
# coldata: metadata table

# Keep only samples that exist in intensity_matrix
coldata <- meta_df %>%
  filter(label %in% colnames(combined_wide)[ecols])

common <- intersect(colnames(combined_wide)[ecols], coldata$label)

coldata <- coldata[match(common, coldata$label), ]



# # Optional: reorder rows to match the column order in intensity_matrix
# coldata <- coldata[match(colnames(combined_wide)[ecols], coldata$label), ]

# # check that sample names match intensity matrix
all(coldata$label %in% colnames(combined_wide)[ecols])

# ----------------------------------------
# Running DEP2
# ----------------------------------------
library(DEP2)

combined_wide$ID <- combined_wide$name
ecols <- grep("^SaMu", colnames(combined_wide))

# Create SummarizedExperiment
coldata$condition <- factor(coldata$sarc_status_bin, levels = c("NoSarc", "Sarc"))
coldata$replicate <- 1

combined_wide <- combined_wide %>% as.data.frame()
rownames(combined_wide) <- combined_wide$ID
combined_wide$name <- combined_wide$ID
coldata <- coldata[, c("label", "condition", "replicate")]
coldata$condition <- as.character(coldata$condition)
# Replace NA with a placeholder like "Unknown"
coldata$condition[is.na(coldata$condition)] <- "Unknown"
coldata$replicate <- as.numeric(gsub("SaMu", "", coldata$label))
se <- make_se(
  combined_wide,
  columns = ecols,
  expdesign = coldata
)

# drop samples with unknown diagnosis
# se <- se[, colData(se)$condition != "Unknown"]
# se <- se[rowData(se)$gene != "NA", ]

# Test different filtering thresholds
# fractions <- seq(0.1, 0.9, by = 0.05)
# thresholds <- c(10, 25, 50, 100)
# 
# param_grid <- expand.grid(
#   fraction = fractions,
#   thr = thresholds
# )

# results <- param_grid %>%
#   rowwise() %>%
#   mutate(
#     n_rows = {
#       se_tmp <- filter_se(
#         se,
#         thr = thr,
#         fraction = fraction
#       )
#       nrow(se_tmp)
#     }
#   ) %>%
#   ungroup()

# p <- ggplot(results, aes(x = fraction, y = n_rows, color = factor(thr))) +
#   geom_line(size = 1) +
#   labs(
#     x = "thr (allowed missing replicates)",
#     y = "missing fraction threshold",
#     color = "missing replicate threshold"
#   ) +
#   theme_minimal()
# ggsave(file.path(parent_dir, "results", "dep2_se_filtering_thresholds.pdf"), p, width = 8, height = 6)

# Filter low-quality proteins
se_filt <- filter_se(se,                       
#                      thr = 50,  ## the threshold of missing number in at least one condition
                     fraction = 0.25 ## the threshold of missing occupancy in each protein
                    )
# NOTE: the documentation for the fraction seems to be incorrect.
# It is not the fraction of missing values, but rather the fraction of valid values required.
print("Dimensions before and after filtering:")
dim(se)
dim(se_filt)

# pdf(file.path(parent_dir, "results", "dep2_filtering_effect.pdf"), width = 20, height = 6)
# plot_frequency(se) + ggtitle("Identification overlap before filter")
# plot_frequency(se_filt) + ggtitle("Identification overlap after before filter")
# dev.off()





# Normalize
n_features <- nrow(assay(se_filt))

if (n_features >= 50) {
  se_norm <- normalize_vsn(se_filt)
} else {
  message("Too few features for VSN (", n_features, "). Using log2 transform instead.")
  assay(se_filt) <- log2(assay(se_filt) + 1)
  se_norm <- se_filt
}
# pdf(file.path(parent_dir, "results", "dep2_normalization_effect.pdf"), width = 6, height = 20)
# plot_normalization(se_filt, se_norm)
# dev.off()



# Plot mean-variance relationship: confirm there is no relationship (horizontal fit)
library(vsn)
# pdf(file.path(parent_dir, "results", "meanSdPlot.pdf"), width = 8, height = 6)
# meanSdPlot(assay(se_norm))
# dev.off()

# Impute missing values
# pdf(file.path(parent_dir, "results", "dep2_missing_values_heatmap.pdf"), width = 8, height = 6)
# plot_missval(se_filt)
# dev.off()

# se_imp <- DEP2::impute(se_norm) # fails
se_imp <- DEP2::impute(se_norm, fun = "MinDet")   # for minimal detection limit
# se_imp <- DEP2::impute(se_norm, fun = "MinProb", q = 0.01)

# Write the imputed, VSN-normalized data to a CSV
# Extract imputed VSN matrix
imputed_matrix <- assay(se_imp)

# transpose so rows = samples, columns = proteins
proteomics_matrix <- t(imputed_matrix)

proteomics_df <- as.data.frame(proteomics_matrix)
proteomics_df$Sample <- rownames(proteomics_df)

# move sample column first
proteomics_df <- proteomics_df[, c("Sample", setdiff(colnames(proteomics_df), "Sample"))]

normed_not_imputed <- assay(se_norm)
write.csv(
  normed_not_imputed,
  file.path(output_dir, paste0("proteomics_", experiment_name, "_vsn_normed_matrix.csv")),
  row.names = TRUE
)

write.csv(
  proteomics_df,
  file.path(output_dir, paste0("proteomics_", experiment_name, "_vsn_imputed_matrix.csv")),
  row.names = FALSE
)

metadata_df <- as.data.frame(colData(se_imp))
write.csv(
  metadata_df,
  file.path(output_dir, paste0("proteomics_", experiment_name, "_vsn_imputed_metadata.csv")),
  row.names = FALSE
)




# -----



# Test differential expression
colData(se_imp)$condition <- factor(colData(se_imp)$condition)
colData(se_imp)$condition <- relevel(
  colData(se_imp)$condition,
  ref = "NoSarc"
)

diff <- test_diff(
  se_imp,
  type = "manual",
  test = "Sarc_vs_NoSarc",
  fdr.type = "BH"
)

# Add significance annotations
results <- add_rejections(diff)
write.csv(as.data.frame(rowData(results)),
          file.path(output_dir, paste0("proteomics_", experiment_name, "_default_results.csv")))

# Quick volcano plot
pdf(file.path(output_dir, paste0("proteomics_", experiment_name, "_default_volcano_plot.pdf")), width = 8, height = 6)
plot_volcano(results, adjusted = F,
             add_threshold_line = "intersect", pCutoff  = 0.05, fcCutoff = 1)
dev.off()

# ----------------------------------------
# Heatmap of the data
# ----------------------------------------
library(ComplexHeatmap)
library(circlize)

# Use normalized intensities
exprs_mat <- assay(se_filt)

log_mat <- log10(exprs_mat + 1e-6)
scaled_mat <- t(scale(t(log_mat)))

# # filter low-variance proteins (for cleaner heatmap)
# # Compute row-wise variance (ignoring NAs)
# num_top_prots <- 100
# row_vars <- apply(exprs_mat, 1, var, na.rm = TRUE)
# # Remove any proteins where variance is NA
# row_vars <- row_vars[!is.na(row_vars)]
# # Get names of top 50 most variable proteins
# top50_proteins <- names(sort(row_vars, decreasing = TRUE))[1:num_top_prots]

# grab by smallest p-value from the DE results
cap_proteins <- FALSE
if (cap_proteins) {
  top50_proteins <- rowData(results) %>%
    as.data.frame() %>%
    arrange(Sarc_vs_NoSarc_p.val) %>%      # sort by this column
    slice_head(n = 50) %>%
    pull(ID)           # extract this column
} else {
  top50_proteins <- rowData(results) %>%
    as.data.frame() %>%
    arrange(Sarc_vs_NoSarc_p.val) %>%      # sort by this column
    pull(ID)           # extract this column
}

# 1. Get the top 50 proteins by smallest p-value
# results_df <- as.data.frame(rowData(results))
# top50_proteins <- results_df %>%
#   arrange(Sarc_vs_NoSarc_p.val) %>%      # or "qval" if that’s the column name
#   slice_head(n = 50) %>%
#   pull(name)              # or "ID" depending on your column

# Subset the expression matrix
# exprs_top50 <- exprs_mat[top50_proteins, ]
scaled_mat_top50 <- scaled_mat[top50_proteins, ]

# Get sample annotation
sample_anno <- data.frame(
  SarcStatus = colData(se_norm)$condition
)
rownames(sample_anno) <- colnames(scaled_mat_top50)

# Compute breakpoints
# Compute breakpoints
sorted_vals <- sort(unique(as.numeric(scaled_mat_top50)), na.last = NA)
min_val <- sorted_vals[1]              # absolute minimum (imputed)
# second_min <- sorted_vals[2]           # second smallest (observed)
max_val <- max(scaled_mat_top50, na.rm = TRUE)  # maximum

# Define color function
col_fun <- colorRamp2(
  breaks = c(min_val, max_val),
  colors = c("black", "red")
)

# scaled_mat_top50[is.na(scaled_mat_top50)] <- 0 # set to 0 for visualization (becomes gray)

# Add sarc/no_sarc annotation
sample_anno <- data.frame(SarcStatus = colData(se_imp)$condition)
rownames(sample_anno) <- colnames(scaled_mat_top50)
ha <- HeatmapAnnotation(
  SarcStatus = sample_anno$SarcStatus,
  col = list(SarcStatus = c("NoSarc" = "steelblue", "Sarc" = "orange", "Unknown" = "gray"))
)

mat_for_clust <- scaled_mat_top50
mat_for_clust[is.na(mat_for_clust)] <- 0

col_order <- order(colnames(scaled_mat_top50))
pdf(file.path( output_dir, paste0("proteomics_", experiment_name, "_heatmap.pdf")), width = 15, height = 12)
Heatmap(
  scaled_mat_top50,
  col = col_fun,
  na_col = "gray80",
  show_row_names = TRUE,                # show protein names
  row_names_side = "left",              # labels on left side
  row_names_gp = gpar(fontsize = 10),   # adjust font size
  name = "scaled Log10(Intensity)",                 # legend title
  top_annotation = ha,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  column_order = col_order,
  # cluster_rows = FALSE,
  # clustering_distance_rows = "euclidean",
  clustering_distance_rows = function(x) dist(mat_for_clust),
  clustering_method_rows = "complete",
  column_title = "Samples",
  # row_title = "Proteins"
)
dev.off()



# Write the cluster assignments to a CSV
peak_cluster_export <- df_long %>%
  group_by(cluster, peak) %>%
  summarise(median_intensity = median(Intensity, na.rm = TRUE), .groups = "drop") %>%
  left_join(peak_clusters_filtered, by = c("cluster", "peak")) %>%
  group_by(cluster) %>%
  mutate(cluster_size = n()) %>%
  ungroup() %>%
  arrange(cluster, desc(median_intensity))

write.csv(
  peak_cluster_export,
  file.path(output_dir, paste0(experiment_name, "_peak_cluster_mapping.csv")),
  row.names = FALSE
)