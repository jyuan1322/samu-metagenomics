# =============================================================================
# 03_lc_feature_exploration.R   (EXPLORATION / one-off, not a routine pipeline)
#
# LC-MS feature exploration for a specific set of input files:
#   1. Audit which features appear in which input files.
#   2. Build a peak-peak correlation network and cluster co-eluting/correlated
#      peaks (igraph connected components above a correlation threshold).
#   3. Collapse each cluster to a representative (mean/max) feature, writing a
#      clustered sample x feature matrix.
#
# The clustered matrix is written in the "wide" layout consumed by
# 02_proteomics_dep2.R (loader = "wide"), so differential expression is run
# there rather than duplicating the DEP2 backbone here.
#
# This script is intentionally kept separate and lightly generalized: it is
# tied to specific inputs (LC_EXPLORE$files) and exploratory choices, not a
# turnkey pipeline. Config in config.R under LC_EXPLORE.
#
# Run:  Rscript 03_lc_feature_exploration.R
# =============================================================================
source("config.R")
source(METADATA_UTILS)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(tidyr)
  library(tibble); library(pheatmap); library(igraph)
})

cfg <- LC_EXPLORE
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Feature presence across files
# ---------------------------------------------------------------------------
feature_file_map <- map_dfr(cfg$files, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  tibble(Metabolite = setdiff(colnames(df), "Sample"), File = basename(file))
})

feature_presence_table <- feature_file_map %>%
  mutate(Present = 1) %>%
  pivot_wider(names_from = File, values_from = Present, values_fill = 0)

print(count(feature_file_map, File, name = "n_features"))
print(feature_presence_table %>%
        mutate(n_files = rowSums(select(., -Metabolite))) %>%
        count(n_files, name = "n_features"))

write_csv(feature_presence_table,
          file.path(cfg$output_dir,
                    paste0(cfg$experiment_name, "_feature_file_presence.csv")))

# ---------------------------------------------------------------------------
# 2. Combine files into one sample x feature matrix (features as rows here)
# ---------------------------------------------------------------------------
combined_wide <- map(cfg$files, function(file) {
  df <- as.data.frame(read_csv(file, show_col_types = FALSE))
  rownames(df) <- df$Sample; df$Sample <- NULL
  df_t <- as.data.frame(t(df)); df_t$name <- rownames(df_t); df_t
}) %>%
  reduce(full_join, by = "name") %>%
  as.data.frame()

colnames(combined_wide) <- gsub("_", "", colnames(combined_wide))  # match labels
rownames(combined_wide) <- combined_wide$name
combined_wide$name <- NULL
combined_wide[is.na(combined_wide)] <- 0

mat <- as.matrix(combined_wide)

# Keep peaks present (non-zero) in at least the configured fraction of samples.
keep <- rowSums(mat != 0) / ncol(mat) >= cfg$peak_filter_threshold
mat <- mat[keep, ]

# ---------------------------------------------------------------------------
# 3. Peak-peak correlation network + clustering
# ---------------------------------------------------------------------------
cor_mat <- cor(t(mat), use = "pairwise.complete.obs")
rownames(cor_mat) <- colnames(cor_mat) <- rownames(mat)

# Correlation heatmap (diagonal blanked) for inspection.
cor_display <- cor_mat; diag(cor_display) <- 0
pdf(file.path(cfg$output_dir, paste0(cfg$experiment_name, "_correlation_heatmap.pdf")),
    width = 20, height = 20)
pheatmap(cor_display, show_rownames = FALSE, show_colnames = FALSE)
dev.off()

# Adjacency: edges where correlation >= threshold.
adj_w <- cor_mat
adj_w[adj_w < cfg$corr_threshold] <- 0
diag(adj_w) <- 0
g <- graph_from_adjacency_matrix(adj_w, mode = "undirected",
                                 weighted = TRUE, diag = FALSE)
comp <- components(g)
peak_clusters <- data.frame(peak = names(comp$membership),
                            cluster = comp$membership, row.names = NULL)

# Drop clusters smaller than the configured minimum size.
sizes <- table(peak_clusters$cluster)
keep_clusters <- names(sizes)[sizes >= cfg$min_cluster_size]
peak_clusters_filtered <- peak_clusters[peak_clusters$cluster %in% keep_clusters, ]

message(sprintf("Samples: %d | filtered peaks: %d | clusters kept: %d | singletons: %d",
                ncol(mat), nrow(mat),
                length(unique(peak_clusters_filtered$cluster)),
                sum(comp$csize == 1)))

# ---------------------------------------------------------------------------
# 4. Collapse each cluster to a representative feature -> wide matrix
# ---------------------------------------------------------------------------
summary_fun <- switch(cfg$cluster_summary,
                      mean = function(x) mean(x, na.rm = TRUE),
                      max  = function(x) max(x, na.rm = TRUE),
                      stop("Unknown LC_EXPLORE$cluster_summary: ", cfg$cluster_summary))

df_clustered <- as.data.frame(mat) %>%
  mutate(peak = rownames(mat)) %>%
  pivot_longer(-peak, names_to = "Sample", values_to = "Intensity") %>%
  inner_join(peak_clusters_filtered, by = "peak") %>%
  group_by(cluster, Sample) %>%
  summarise(Intensity = summary_fun(Intensity), .groups = "drop")

# Output layout: Sample column + one column per cluster feature, matching the
# "wide" loader of 02_proteomics_dep2.R.
df_wide <- df_clustered %>%
  mutate(cluster = paste0("Cluster", cluster)) %>%
  pivot_wider(names_from = cluster, values_from = Intensity)

out_csv <- file.path(cfg$output_dir,
                     paste0(cfg$experiment_name, "_clustered_features.csv"))
write_csv(df_wide, out_csv)
message("Wrote clustered feature matrix: ", out_csv)
message("Feed this into 02_proteomics_dep2.R with loader = \"wide\" for DE.")
message("03_lc_feature_exploration.R complete.")
