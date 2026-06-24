# =============================================================================
# 01_nmr_quorum.R
#
# NMR metabolite and quorum-sensing-molecule analysis: join assay tables to
# metadata, log-transform, scaled heatmaps, PCA/PCoA, and per-metabolite
# Wilcoxon tests with boxplots. Config in config.R under NMR.
#
# Run:  Rscript 01_nmr_quorum.R
# =============================================================================
source("config.R")
source(METADATA_UTILS)

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr); library(ggplot2)
  library(pheatmap); library(tibble)
  library(phyloseq); library(vegan)
})

cfg <- NMR
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(cfg$output_dir)

# ---------------------------------------------------------------------------
# Load metadata + assay tables, build a merged frame
# ---------------------------------------------------------------------------
# Metadata is recoded from the raw CSV (not the metagenomics output) so that
# NMR/quorum samples outside the metagenomics cohort are still covered.
meta_df <- load_samu_metadata(
  METADATA_CSV,
  keep_cols    = META_KEEP_COLS,
  fullsamu_col = FULLSAMU_COL,
  group_levels = GROUP_LEVELS,
  recode_smoke_alcohol = RECODE_SMOKE_ALCOHOL
)

nmr_concs    <- read.csv(file.path(cfg$input_dir, cfg$nmr_csv))
quorum_concs <- read.csv(file.path(cfg$input_dir, cfg$quorum_csv))

nmr_cols    <- colnames(nmr_concs)[cfg$nmr_first_col:ncol(nmr_concs)]
quorum_cols <- str_subset(colnames(quorum_concs), cfg$quorum_col_regex)

# Below-detection-limit placeholder -> 0 in NMR columns.
nmr_concs[nmr_cols] <- lapply(nmr_concs[nmr_cols], function(x) {
  x[x == cfg$nmr_below_detection] <- "0"; x
})

# load_samu_metadata already applied the Full.SaMu filter; a left join keeps
# every retained metadata sample and attaches assay values where present.
merged_df <- meta_df %>%
  merge(nmr_concs[, c("record_id", nmr_cols)],   by = "record_id", all.x = TRUE) %>%
  merge(quorum_concs[, c("record_id", quorum_cols)], by = "record_id", all.x = TRUE)

merged_df[[GROUP_VAR]] <- factor(merged_df[[GROUP_VAR]], levels = GROUP_LEVELS)

write.csv(merged_df, "nmr_quorum_FullSaMu.csv", row.names = FALSE, quote = TRUE)

# ---------------------------------------------------------------------------
# Helper: build a numeric sample x feature matrix from selected columns,
# dropping all-NA samples. record_id becomes the row name.
# ---------------------------------------------------------------------------
build_matrix <- function(df, cols) {
  m <- df %>%
    mutate(across(all_of(cols), as.numeric)) %>%  # NAs-by-coercion expected
    select(all_of(cols)) %>%
    as.matrix()
  rownames(m) <- df$record_id
  m[rowSums(is.na(m)) < ncol(m), , drop = FALSE]
}

# Helper: ordered group annotation + the row gap between groups, aligned to a
# matrix's row order.
group_annotation <- function(mat, df) {
  anno <- df %>%
    filter(record_id %in% rownames(mat)) %>%
    arrange(match(record_id, rownames(mat))) %>%
    select(record_id, all_of(GROUP_VAR)) %>%
    column_to_rownames("record_id")
  anno[[GROUP_VAR]] <- factor(anno[[GROUP_VAR]], levels = GROUP_LEVELS)
  ord <- order(anno[[GROUP_VAR]])
  list(anno = anno[ord, , drop = FALSE], order = ord,
       gap = sum(anno[[GROUP_VAR]] == GROUP_LEVELS[1]))
}

heat_colors <- colorRampPalette(c("black", "red"))(100)

# ---------------------------------------------------------------------------
# NMR: log10 concentrations (NOTE: NMR is not compositional, so we use log
# concentrations directly rather than relative abundances).
# ---------------------------------------------------------------------------
nmr_mat <- build_matrix(merged_df, nmr_cols)
log_nmr <- log10(nmr_mat + cfg$epsilon)

ga <- group_annotation(log_nmr, merged_df)
log_nmr_scaled <- scale(log_nmr, center = TRUE, scale = TRUE)[ga$order, , drop = FALSE]

write.csv(log_nmr_scaled, "samu_nmr_logrel_scaled.csv", row.names = TRUE, quote = TRUE)

pheatmap(log_nmr_scaled, color = heat_colors, scale = "none",
         clustering_method = "ward.D2", cluster_cols = TRUE, cluster_rows = FALSE,
         annotation_row = ga$anno, gaps_row = ga$gap, fontsize = 8,
         border_color = NA, filename = "samu_nmr_logrel_scaled.pdf",
         width = 12, height = 12)

# ---------------------------------------------------------------------------
# NMR beta diversity (Bray-Curtis). NOTE: Bray-Curtis is only a rough choice
# here since NMR data are not compositional; interpret with caution.
# ---------------------------------------------------------------------------
otu_mat <- otu_table(nmr_mat, taxa_are_rows = FALSE)
ps_metab <- phyloseq(otu_mat, sample_data(ga$anno))
ps_metab <- prune_taxa(!(taxa_names(ps_metab) %in% cfg$beta_remove_metabs),
                       ps_metab)
m <- as(otu_table(ps_metab), "matrix"); m[is.na(m)] <- 0
otu_table(ps_metab) <- otu_table(m, taxa_are_rows = taxa_are_rows(ps_metab))

bc_dist  <- phyloseq::distance(ps_metab, method = "bray")
pcoa_res <- ordinate(ps_metab, method = "PCoA", distance = bc_dist)

p <- plot_ordination(ps_metab, pcoa_res, color = GROUP_VAR, shape = "sex") +
  geom_point(size = 3) + theme_classic() +
  labs(title = "PCoA of Metabolite Bray-Curtis Distances")
ggsave("samu_nmr_beta_diversity_pcoa.pdf", p, width = 6, height = 5)

# PCoA colored by total concentration ("read depth" analogue).
read_depth <- sample_sums(ps_metab)
p2 <- plot_ordination(ps_metab, pcoa_res) +
  geom_point(aes(color = read_depth[rownames(pcoa_res$vectors)]), size = 3) +
  scale_color_viridis_c() + theme_classic() +
  labs(title = "PCoA Colored by Sum of Concentrations (microM)", color = "Total conc")
ggsave("samu_nmr_beta_diversity_pcoa_color_by_conc_depth.pdf", p2, width = 6, height = 5)

# ---------------------------------------------------------------------------
# NMR PCA (scaled log concentrations), colored by group / age / depth
# ---------------------------------------------------------------------------
log_nmr_scaled_clean <- log_nmr_scaled[, colSums(is.na(log_nmr_scaled)) == 0,
                                       drop = FALSE]
pca_res <- prcomp(log_nmr_scaled_clean, center = FALSE, scale. = FALSE)
pca_df  <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- cbind(pca_df,
                merged_df[match(pca_df$Sample, merged_df$record_id),
                          c(GROUP_VAR, "age_def"), drop = FALSE])
pca_df$read_depth <- read_depth[pca_df$Sample]

ggsave("samu_nmr_pca_1_color_by_sarc_status.pdf",
       ggplot(pca_df, aes(PC1, PC2, color = .data[[GROUP_VAR]])) +
         geom_point(size = 3) + theme_classic() +
         labs(title = "PCA colored by group", color = "Group"),
       width = 6, height = 5)
ggsave("samu_nmr_pca_color_by_age.pdf",
       ggplot(pca_df, aes(PC1, PC2, color = age_def)) +
         geom_point(size = 3) + scale_color_viridis_c() + theme_classic() +
         labs(title = "PCA colored by age", color = "Age (years)"),
       width = 6, height = 5)
ggsave("samu_nmr_pca_color_by_conc_depth.pdf",
       ggplot(pca_df, aes(PC1, PC2, color = read_depth)) +
         geom_point(size = 3) + scale_color_viridis_c(option = "C") +
         theme_classic() +
         labs(title = "PCA colored by total concentration", color = "Total conc"),
       width = 6, height = 5)

# ---------------------------------------------------------------------------
# Quorum sensing: compositional, so use log relative abundances
# ---------------------------------------------------------------------------
quorum_mat <- build_matrix(merged_df, quorum_cols)
quorum_mat[is.na(quorum_mat)] <- 0                       # unique to quorum data
quorum_mat <- quorum_mat[rowSums(quorum_mat) > 0, , drop = FALSE]
quorum_mat <- quorum_mat[, colSums(quorum_mat) > 0, drop = FALSE]

relab <- sweep(quorum_mat, 1, rowSums(quorum_mat), "/")
log_quorum <- log10(relab + cfg$epsilon)

gaq <- group_annotation(log_quorum, merged_df)
log_quorum <- log_quorum[gaq$order, , drop = FALSE]

# Gray for the minimum (below-detection) bin, red ramp above.
q_colors <- c("gray80", colorRampPalette(c("black", "red"))(99))
q_breaks <- seq(min(log_quorum, na.rm = TRUE),
                max(log_quorum, na.rm = TRUE) + 1e-12, length.out = length(q_colors))
pheatmap(log_quorum, color = q_colors, breaks = q_breaks, scale = "none",
         clustering_method = "ward.D2", cluster_cols = TRUE, cluster_rows = TRUE,
         annotation_row = gaq$anno, fontsize = 8, border_color = NA,
         filename = "samu_quorum_logrel_unscaled.pdf", width = 12, height = 12)

# ---------------------------------------------------------------------------
# Per-metabolite Wilcoxon tests + boxplots (NMR metabolites)
# ---------------------------------------------------------------------------
long_df <- merged_df %>%
  select(all_of(GROUP_VAR), all_of(nmr_cols)) %>%
  pivot_longer(all_of(nmr_cols), names_to = "metabolite", values_to = "value") %>%
  mutate(value = as.numeric(value))

pvals <- long_df %>%
  group_by(metabolite) %>%
  summarise(p_value = wilcox.test(value ~ .data[[GROUP_VAR]])$p.value,
            .groups = "drop") %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

y_positions <- long_df %>%
  mutate(log_value = log10(value + 0.01)) %>%
  group_by(metabolite) %>%
  summarise(y_pos = max(log_value, na.rm = TRUE) * 1.1, .groups = "drop")

pvals_plot <- pvals %>%
  left_join(y_positions, by = "metabolite") %>%
  mutate(label = paste0("p = ", signif(p_value, 2),
                        "\n(padj = ", signif(p_adj, 2), ")"))

# Order metabolites by raw p-value on the x-axis.
lvl <- pvals_plot %>% arrange(p_value) %>% pull(metabolite)
long_df$metabolite    <- factor(long_df$metabolite, levels = lvl)
pvals_plot$metabolite <- factor(pvals_plot$metabolite, levels = lvl)

p <- ggplot(long_df, aes(metabolite, log10(value + 1), fill = .data[[GROUP_VAR]])) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
  geom_jitter(aes(color = .data[[GROUP_VAR]]),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.5, size = 0.7, show.legend = FALSE) +
  geom_text(data = pvals_plot, aes(metabolite, y_pos, label = label),
            inherit.aes = FALSE, size = 3.5) +
  theme_classic() +
  labs(x = "Metabolite", y = "log10(Concentration + 1)", fill = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("nmr_wilcoxon_boxplot.pdf", p, width = 30, height = 6, dpi = 300)

message("01_nmr_quorum.R complete.")
