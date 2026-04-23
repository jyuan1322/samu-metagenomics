library(dplyr)
library(stringr)
library(pheatmap)

input_dir <- "/data/local/jy1008/SaMu"
meta_df <- read.csv(file.path(input_dir, "results", "02102026", "meta_df_FullSaMu.csv"))
read_counts <- read.csv(file.path(input_dir, "results", "02102026", "meta_filtered.csv"))
nmr_concs <- read.csv(file.path(input_dir, "metadata", "20251027_SaMu_NMR_ureum_creat_v10.csv"))
quorum_concs <- read.csv(file.path(input_dir, "metadata", "20251030_SaMu_QSP_data_v1.csv"))

nmr_cols <- colnames(nmr_concs)[12:ncol(nmr_concs)]
quorum_cols <- colnames(quorum_concs) %>%
  str_subset("_Quant\\.Prob$")
meta_cols <- c("record_id", "Block", "sarc_status_bin", "age_def", "sex", "smke", "alco", "nutr_score", "bmi", "stvol")

# Replace "UNK_E_1" with 0 (below detection limit) in NMR data
nmr_concs[nmr_cols] <- lapply(nmr_concs[nmr_cols], function(x) {
  x[x == "UNK_E_1"] <- "0"
  x
})

merged_df <- merge(
  x = meta_df,
  y = nmr_concs[, c("record_id", nmr_cols)],
  by = "record_id",
  all = TRUE # or TRUE, depending on the join type
)
merged_df <- merge(
  x = merged_df,
  y = quorum_concs[, c("record_id", quorum_cols)],
  by = "record_id",
  all = TRUE # or TRUE, depending on the join type
)

# filter out not full SaMu samples
merged_df <- merged_df %>%
  filter(
    !is.na(Full.SaMu) &
    Full.SaMu == 1
  ) %>% as.data.frame()

# write matrix
write.csv(
  merged_df,
  "nmr_quorum_FullSaMu.csv",
  row.names = FALSE,
  quote = TRUE
)

# create a heatmap of the NMR concentrations
# 1. Convert selected NMR columns to numeric and turn into matrix
metabolite_matrix <- merged_df %>%
  mutate(across(all_of(nmr_cols), as.numeric)) %>%  # safe numeric conversion
  select(all_of(nmr_cols)) %>%
  as.matrix()

# Add sample IDs as rownames
rownames(metabolite_matrix) <- merged_df$record_id

# 2. Remove rows that are all NA
metabolite_matrix <- metabolite_matrix[
  rowSums(is.na(metabolite_matrix)) < ncol(metabolite_matrix),
  ,
  drop = FALSE
]
# This warning is expected if there are non-numeric entries in the data
# Caused by warning:
# ! NAs introduced by coercion
# ℹ Run `dplyr::last_dplyr_warnings()` to see the 39 remaining warnings.

# NOTE (2/11/2026): For NMR data, do not use relative abundance, only log concs
# # 3. Compute per-sample sums
# sum_per_sample <- rowSums(metabolite_matrix, na.rm = TRUE)

# # 4. Compute relative abundances
# relabund_matrix <- sweep(
#   metabolite_matrix,
#   MARGIN = 1,
#   STATS = sum_per_sample,
#   FUN = "/"
# )

# 5. Log-transform with small epsilon
epsilon <- 1e-6
# logrel_matrix <- log10(relabund_matrix + epsilon)
logrel_matrix <- log10(metabolite_matrix + epsilon)

metadata_df_cleaned <- merged_df %>%
  filter(record_id %in% rownames(logrel_matrix)) %>%
  arrange(match(record_id, rownames(logrel_matrix)))
# verify matching entries
all(metadata_df_cleaned$record_id == rownames(logrel_matrix))


# annotation_df <- metadata_df_cleaned %>%
#   select(record_id, Block, sarc_status_bin, sex, smke, age_def, nutr_score, bmi) %>%   # <--- choose your columns
#   tibble::column_to_rownames("record_id")
# # ensure rownames match
# annotation_df <- annotation_df[rownames(logrel_matrix), , drop = FALSE]
# annotation_df$Block <- as.factor(annotation_df$Block)
# annotation_df$sarc_status_bin <- factor(annotation_df$sarc_status_bin, levels = c("NoSarc", "Sarc"))
# annotation_df$sex <- as.factor(annotation_df$sex)
# annotation_df$smke <- as.factor(annotation_df$smke)
# plot only sarc_status_bin
annotation_df <- metadata_df_cleaned %>%
  select(record_id, sarc_status_bin) %>%
  tibble::column_to_rownames("record_id") %>%
  as.data.frame()

annotation_df$sarc_status_bin <- factor(annotation_df$sarc_status_bin,
                                         levels = c("NoSarc", "Sarc"))

# Z-score each metabolite (column)
logrel_matrix_scaled <- scale(logrel_matrix, center = TRUE, scale = TRUE)

# Keep rownames and column names
rownames(logrel_matrix_scaled) <- rownames(logrel_matrix)
colnames(logrel_matrix_scaled) <- colnames(logrel_matrix)

# Order samples by sarc status
sample_order <- order(annotation_df$sarc_status_bin)

logrel_matrix_scaled <- logrel_matrix_scaled[sample_order, , drop = FALSE]
annotation_df <- annotation_df[sample_order, , drop = FALSE]

gaps_row <- sum(annotation_df$sarc_status_bin == "NoSarc")

write.csv(
  logrel_matrix_scaled,
  "samu_nmr_logrel_scaled.csv",
  row.names = TRUE,
  quote = TRUE
)

my_colors <- colorRampPalette(c("black", "red"))(100)
pheatmap(
  logrel_matrix_scaled,
  color = my_colors,
  scale = "none",
  clustering_method = "ward.D2",
  cluster_cols = TRUE,
  cluster_rows = FALSE,
  annotation_row = annotation_df,   # <-- sample annotations
  gaps_row = gaps_row,
  fontsize = 8,
  border_color = NA,
  filename = "samu_nmr_logrel_scaled.pdf",
  width = 12,
  height = 12
)

# Z-score each metabolite (column)
metabolite_matrix_scaled <- scale(metabolite_matrix, center = TRUE, scale = TRUE)

# Keep rownames and column names
rownames(metabolite_matrix_scaled) <- rownames(metabolite_matrix_scaled)
colnames(metabolite_matrix_scaled) <- colnames(metabolite_matrix_scaled)

pheatmap(
  metabolite_matrix_scaled,
  color = my_colors,
  scale = "none",
  clustering_method = "ward.D2",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  annotation_row = annotation_df,   # <-- sample annotations
  fontsize = 8,
  border_color = NA,
  filename = "samu_nmr_concs_scaled.pdf",
  width = 12,
  height = 12
)


# NOTE (2/11/2026): Bray-Curtis is not appropriate for NMR because it's not compositional
# beta diversity of NMR concentrations
library(phyloseq)
library(vegan)
library(ggplot2)

# Convert metabolite matrix to an OTU table
# taxa_are_rows = FALSE because rows = samples, cols = metabolites
otu_mat <- otu_table(metabolite_matrix, taxa_are_rows = FALSE)

# Convert annotation to sample_data
sample_data_obj <- sample_data(annotation_df)

# Create phyloseq object
ps_metab <- phyloseq(otu_mat, sample_data_obj)

# Remove samples with any NA values in the OTU table
# ps_metab_rel <- prune_samples(
#   rowSums(is.na(otu_table(ps_metab_rel))) == 0,
#   ps_metab_rel
# )
# Remove metabolites with any NA values across metabolites
# ps_metab <- prune_taxa(
#   colSums(is.na(otu_table(ps_metab))) == 0,
#   ps_metab
# )
# Remove particular metabolites columns with high missingness
# remove_metabs <- c("Glucose", "Galactose", "Urea", "Cadaverine")
remove_metabs <- c("Glucose", "Galactose")
ps_metab <- prune_taxa(
  !(taxa_names(ps_metab) %in% remove_metabs),
  ps_metab
)
# Replace any remaining NAs with zeros
# Extract OTU table as a matrix
otu_mat <- as(otu_table(ps_metab), "matrix")

# Replace NA with 0
otu_mat[is.na(otu_mat)] <- 0

# Put it back into the phyloseq object
otu_table(ps_metab) <- otu_table(
  otu_mat,
  taxa_are_rows = taxa_are_rows(ps_metab)
)
any(is.na(otu_table(ps_metab)))
# should be FALSE



# optional: convert to relative abundance
# ps_metab_rel <- transform_sample_counts(ps_metab, function(x) x / sum(x))

bc_dist <- phyloseq::distance(ps_metab, method = "bray")

pcoa_res <- ordinate(ps_metab, method = "PCoA", distance = bc_dist)

# Use phyloseq's built-in plotting
p <- plot_ordination(ps_metab, pcoa_res, color = "sarc_status_bin", shape = "sex") +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "PCoA of Metabolite Bray-Curtis Distances")
ggsave("samu_nmr_beta_diversity_pcoa.pdf", plot = p, width = 6, height = 5)

# check if read depth is responsible for PC1
pcoa_df <- cbind(
  sample_data(ps_metab) %>% as.data.frame(),
  pcoa_res$vectors[, 1:2]    # PC1 and PC2
)
# Name the columns
colnames(pcoa_df)[(ncol(pcoa_df)-1):ncol(pcoa_df)] <- c("PC1", "PC2")

read_depth <- sample_sums(ps_metab)
pcoa_df$Sample <- rownames(pcoa_df)
pcoa_df$read_depth <- read_depth[pcoa_df$Sample]

p2 <- plot_ordination(ps_metab, pcoa_res) +
  geom_point(aes(color = read_depth), size = 3) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(title = "PCoA Colored by Sum of Concentrations in microM")
ggsave("samu_nmr_beta_diversity_pcoa_color_by_conc_depth.pdf", plot = p2, width = 6, height = 5)


# PCA of scaled log relative abundances
logrel_matrix_scaled_clean <- logrel_matrix_scaled[, colSums(is.na(logrel_matrix_scaled)) == 0, drop = FALSE]
pca_res <- prcomp(logrel_matrix_scaled_clean, center = FALSE, scale. = FALSE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- cbind(pca_df, annotation_df[rownames(pca_df), , drop = FALSE])

read_depth <- sample_sums(ps_metab)  # sums of counts per sample
pca_df$read_depth <- read_depth[rownames(pca_df)]


p1 <- ggplot(pca_df, aes(PC1, PC2, color = sarc_status_bin)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "PCA colored by sarcopenia status", color = "Sarcopenia Status")
ggsave("samu_nmr_pca_1_color_by_sarc_status.pdf", plot = p1, width = 6, height = 5)

p2 <- ggplot(pca_df, aes(PC1, PC2, color = age_def)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(title = "PCA colored by age", color = "Age (years)")
ggsave("samu_nmr_pca_color_by_age.pdf", plot = p2, width = 6, height = 5)

p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = read_depth)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "C") +  # continuous color scale
  theme_classic() +
  labs(title = "PCA colored by scaled log concentration (microM)", color = "Read Depth")
ggsave("samu_nmr_pca_color_by_conc_depth.pdf", plot = p3, width = 6, height = 5)





# create a heatmap of the Quorum sensing molecule concentrations






# create a heatmap of the quorum sensing measurements
# 1. Convert selected quorum columns to numeric and turn into matrix
quorum_matrix <- merged_df %>%
  mutate(across(all_of(quorum_cols), as.numeric)) %>%  # safe numeric conversion
  select(all_of(quorum_cols)) %>%
  as.matrix()

# Add sample IDs as rownames
rownames(quorum_matrix) <- merged_df$record_id

# 2. Remove rows that are all NA
quorum_matrix <- quorum_matrix[
  rowSums(is.na(quorum_matrix)) < ncol(quorum_matrix),
  ,
  drop = FALSE
]
# This warning is expected if there are non-numeric entries in the data
# Caused by warning:
# ! NAs introduced by coercion
# ℹ Run `dplyr::last_dplyr_warnings()` to see the 39 remaining warnings.

# 3. Convert remaining NA values to 0 (unique to quorum data)
quorum_matrix[is.na(quorum_matrix)] <- 0

# 3. Compute per-sample sums
sum_per_sample <- rowSums(quorum_matrix, na.rm = TRUE)
# remove samples with zero sum to avoid division by zero
quorum_matrix <- quorum_matrix[sum_per_sample > 0, , drop = FALSE]
# remove columns with zero sum to avoid division by zero
col_sums <- colSums(quorum_matrix, na.rm = TRUE)
quorum_matrix <- quorum_matrix[, col_sums > 0, drop = FALSE]

# 4. Compute relative abundances
sum_per_sample <- rowSums(quorum_matrix, na.rm = TRUE)
relabund_matrix <- sweep(
  quorum_matrix,
  MARGIN = 1,
  STATS = sum_per_sample,
  FUN = "/"
)

# 5. Log-transform with small epsilon
epsilon <- 1e-6
logrel_matrix <- log10(relabund_matrix + epsilon)

metadata_df_cleaned <- merged_df %>%
  filter(record_id %in% rownames(logrel_matrix)) %>%
  arrange(match(record_id, rownames(logrel_matrix)))
# verify matching entries
all(metadata_df_cleaned$record_id == rownames(logrel_matrix))


# annotation_df <- metadata_df_cleaned %>%
#   select(record_id, Block, sarc_status_bin, sex, smke, age_def, nutr_score, bmi) %>%   # <--- choose your columns
#   tibble::column_to_rownames("record_id")
# # ensure rownames match
# annotation_df <- annotation_df[rownames(logrel_matrix), , drop = FALSE]
# annotation_df$Block <- as.factor(annotation_df$Block)
# annotation_df$sarc_status_bin <- factor(annotation_df$sarc_status_bin, levels = c("NoSarc", "Sarc"))
# annotation_df$sex <- as.factor(annotation_df$sex)
# annotation_df$smke <- as.factor(annotation_df$smke)

annotation_df <- metadata_df_cleaned %>%
  select(record_id, sarc_status_bin) %>%
  tibble::column_to_rownames("record_id") %>%
  as.data.frame()

annotation_df$sarc_status_bin <- factor(annotation_df$sarc_status_bin,
                                         levels = c("NoSarc", "Sarc"))

sample_order <- order(annotation_df$sarc_status_bin)

logrel_matrix <- logrel_matrix[sample_order, , drop = FALSE]
annotation_df <- annotation_df[sample_order, , drop = FALSE]

gaps_row <- sum(annotation_df$sarc_status_bin == "NoSarc")

# Z-score each metabolite (column)
# logrel_matrix_scaled <- scale(logrel_matrix, center = TRUE, scale = TRUE)

# Keep rownames and column names
# rownames(logrel_matrix_scaled) <- rownames(logrel_matrix)
# colnames(logrel_matrix_scaled) <- colnames(logrel_matrix)

# set the minimum to be gray
# my_colors <- colorRampPalette(c("black", "red"))(100)
# 1. Define original palette (excluding gray)
orig_colors <- colorRampPalette(c("black", "red"))(99)
# 2. Prepend gray for minimum
my_colors <- c("gray80", orig_colors)
# 3. Get min/max of the matrix (or log matrix)
min_val <- min(logrel_matrix, na.rm = TRUE)
max_val <- max(logrel_matrix, na.rm = TRUE)
# 4. Create strictly increasing breaks
# length(breaks) = length(colors)
my_breaks <- seq(min_val, max_val + 1e-12, length.out = length(my_colors))
pheatmap(
  # logrel_matrix_scaled,
  logrel_matrix,
  color = my_colors,
  breaks = my_breaks,
  scale = "none",
  clustering_method = "ward.D2",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  annotation_row = annotation_df,   # <-- sample annotations
  fontsize = 8,
  border_color = NA,
  filename = "samu_quorum_logrel_unscaled.pdf",
  width = 12,
  height = 12
)

# Z-score each metabolite (column)
quorum_matrix_scaled <- scale(quorum_matrix, center = TRUE, scale = TRUE)

# Keep rownames and column names
rownames(quorum_matrix_scaled) <- rownames(quorum_matrix_scaled)
colnames(quorum_matrix_scaled) <- colnames(quorum_matrix_scaled)

pheatmap(
  quorum_matrix_scaled,
  color = my_colors,
  scale = "none",
  clustering_method = "ward.D2",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  annotation_row = annotation_df,   # <-- sample annotations
  fontsize = 8,
  border_color = NA,
  # filename = "samu_nmr_concs_scaled.pdf",
  width = 12,
  height = 12
)





# 2/17/2026: Wilcoxon rank-sum test and boxplots for each metabolite
library(dplyr)
library(tidyr)
library(ggplot2)

merged_df$sarc_status_bin <- factor(merged_df$sarc_status_bin, levels = c("NoSarc", "Sarc"))

long_df <- merged_df %>%
  select(sarc_status_bin, all_of(nmr_cols)) %>%
  pivot_longer(
    cols = all_of(nmr_cols),
    names_to = "metabolite",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(value))


summary_df <- long_df %>%
  group_by(metabolite, sarc_status_bin) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# order metabolites by absolute difference in median between sarc and no sarc
order_df <- summary_df %>%
  select(metabolite, sarc_status_bin, median) %>%
  pivot_wider(names_from = sarc_status_bin, values_from = median) %>%
  mutate(diff = abs(Sarc - NoSarc)) %>%
  arrange(desc(diff))


summary_df$metabolite <- factor(
  summary_df$metabolite,
  levels = order_df$metabolite
)

pvals <- long_df %>%
  group_by(metabolite) %>%
  summarise(
    p_value = wilcox.test(value ~ sarc_status_bin)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

summary_df <- summary_df %>%
  left_join(pvals, by = "metabolite")




# boxplot with p-values
pvals_plot <- summary_df %>%
  select(metabolite, p_value, p_adj) %>%
  distinct()

y_positions <- long_df %>%
  mutate(log_value = log10(value + 0.01)) %>%
  group_by(metabolite) %>%
  summarise(
    y_pos = max(log_value, na.rm = TRUE) * 1.1,
    .groups = "drop"
  )

pvals_plot <- left_join(pvals_plot, y_positions, by = "metabolite")
pvals_plot <- pvals_plot %>%
  mutate(label = paste0("p = ", signif(p_value, 2), "\n(padj = ", signif(p_adj, 2), ")"))

# reorder x-axis labels according to p-value
long_df$metabolite <- factor(
  long_df$metabolite,
  levels = pvals_plot %>% arrange(p_value) %>% pull(metabolite)
)
pvals_plot$metabolite <- factor(
  pvals_plot$metabolite,
  levels = levels(long_df$metabolite)
)

p <- ggplot(long_df,
       aes(x = metabolite,
           y = log10(value + 1),
           fill = sarc_status_bin)) +
  geom_boxplot(position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_jitter(aes(color = sarc_status_bin),
              position = position_jitterdodge(jitter.width = 0.2,
                                              dodge.width = 0.8),
              alpha = 0.5,
              size = 0.7,
              show.legend = FALSE) +
  geom_text(data = pvals_plot,
            aes(x = metabolite,
                y = y_pos,
                label = label),
            inherit.aes = FALSE,
            size = 3.5) +
  theme_classic() +
  labs(
    x = "Metabolite",
    y = "log10(Concentration + 0.01)",
    fill = "Sarc Status"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("nmr_wilcoxon_boxplot_02172026.pdf",
       plot = p,
       width = 30,
       height = 6,
       dpi = 300)