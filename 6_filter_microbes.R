library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(tools)
library(stringr)
library(forcats)
library(viridis)
library(vegan)
library(phyloseq)
library(DESeq2)

input_dir <- "/data/local/jy1008/SaMu/metaphlan_out_09172025/all_merged_fastqs"

# metadata
meta_df <- read.csv("/data/local/jy1008/SaMu/SaMu_sarcopeniestatus_majorcovariates_v11_metadata.csv",
                header = TRUE, stringsAsFactors = FALSE)
# create link to files
meta_df$File_ID <- NA_character_
colname <- "block1.20231011_AV224503_4520_1.RawData.4520.tar"
# find non-empty rows
nonempty_idx <- trimws(meta_df[[colname]]) != ""
# update only those rows
meta_df$File_ID[nonempty_idx] <- stringr::str_extract(meta_df[[colname]][nonempty_idx], "SaMu[0-9]+")
colname <- "block2.20240327_AV224503_4734_1.RawData.4734.tar"
nonempty_idx <- trimws(meta_df[[colname]]) != ""
meta_df$File_ID[nonempty_idx] <- stringr::str_extract(meta_df[[colname]][nonempty_idx], "\\d+_\\d+")
colname <- "block5.20250124_AV242402_4915_1.RawData.4915.tar"
nonempty_idx <- trimws(meta_df[[colname]]) != ""
meta_df$File_ID[nonempty_idx] <- stringr::str_extract(meta_df[[colname]][nonempty_idx], "\\d+_\\d+_\\d+_SaMu\\d+")
meta_df <- meta_df[!is.na(meta_df$File_ID), ]
rownames(meta_df) <- meta_df$File_ID

# List all MetaPhlAn output files
files <- list.files(input_dir, pattern = "_profile.txt$", full.names = TRUE)

# Function to read one MetaPhlAn file
read_metaphlan <- function(file) {
  # Read all lines
  lines <- readLines(file)
  
  # Find the header line (starts with "#clade_name")
  header_line <- lines[grepl("^#clade_name", lines)]
  
  # Remove leading "#", split by whitespace
  header <- strsplit(sub("^#", "", header_line), "\\s+")[[1]]
  
  # Read the table skipping all comment lines
  df <- read_tsv(file, comment = "#", col_names = header, col_types = cols())
  
  # Add sample name
  df$Sample <- file_path_sans_ext(basename(file))
  
  # Keep only the columns we need
  df %>%
    select(Sample, clade_name, relative_abundance, estimated_number_of_reads_from_the_clade)
}


# Combine all samples into one dataframe
comb_df <- do.call(rbind, lapply(files, read_metaphlan))

# Convert relative abundance to numeric
comb_df <- comb_df %>%
  mutate(relative_abundance = as.numeric(relative_abundance))

# Keep only species-level groups
meta_species <- comb_df %>%
  dplyr::filter(str_detect(clade_name, "s__")) %>%  # keep species-level entries
  dplyr::filter(!str_detect(clade_name, "\\|t__")) %>% # drop strains
  dplyr::rename(Species = clade_name) %>%
  group_by(Sample, Species) %>%
  summarise(
    relative_abundance = sum(relative_abundance, na.rm = TRUE),
    estimated_number_of_reads_from_the_clade = sum(estimated_number_of_reads_from_the_clade, na.rm = TRUE),
    .groups = "drop"
  )

# Check that sums of relative abundances are near 100
temp <- meta_species %>%
  group_by(Sample) %>%
  summarise(total_rel_abundance = sum(relative_abundance, na.rm = TRUE))
print("Sums of relative abundances at species level")
print(temp)

# Elbow plot
# Define a sequence of abundance thresholds, e.g., 0% to 1%
abundance_thresholds <- seq(0, 2, by = 0.02)  # 0% to 2%
prevalence_thresholds <- c(1, 5, 10, 20)  # number of samples a species must appear in

# Build elbow data for multiple prevalence thresholds
elbow_df <- expand_grid(
  threshold = abundance_thresholds,
  min_samples = prevalence_thresholds
) %>%
  mutate(
    CladesRetained = map2_int(threshold, min_samples, ~ {
      meta_species %>%
        group_by(Species) %>%
        summarize(n_samples_above = sum(relative_abundance >= .x), .groups = "drop") %>%
        filter(n_samples_above >= .y) %>%
        nrow()
    })
  )

# Plot multiple lines
p <- ggplot(elbow_df, aes(x = threshold, y = CladesRetained, color = factor(min_samples))) +
  geom_line(size = 1.2) +
  scale_x_continuous(
    breaks = seq(0, 2, by = 0.1)
  ) +
  theme_minimal() +
  labs(
    x = "Relative Abundance Threshold (%)",
    y = "Number of Species Passing Filter",
    color = "Min Samples Passing Threshold",
    title = "Elbow Plot for Abundance Filtering"
  )
ggsave("elbow_plot_species_level_09232025.pdf", plot = p, width = 10, height = 6)


# ----------



# Identify species that pass the threshold in at least 20 samples
threshold = 0.1
min_samples = 20
species_keep <- meta_species %>%
  group_by(Species) %>%
  summarize(n_samples_above_threshold = sum(relative_abundance >= threshold)) %>%
  filter(n_samples_above_threshold >= min_samples) %>%
  pull(Species)

# Filter original dataframe
meta_filtered <- meta_species %>%
  filter(Species %in% species_keep)

# Print number of clades before and after filtering
n_before <- length(unique(meta_species$Species))
n_after <- length(unique(meta_filtered$Species))
cat("Number of clades before filtering:", n_before, "\n")
cat("Number of clades after filtering:", n_after, "\n")
# Number of clades before filtering: 2128 
# Number of clades after filtering: 785

# after performing filtering, renormalize remaining species' rel abundances to sum to 1
meta_filtered <- meta_filtered %>%
  group_by(Sample) %>%
  mutate(relative_abundance = relative_abundance / sum(relative_abundance) * 100) %>%
  ungroup()

# Check that sums of relative abundances are near 100
temp <- meta_filtered %>%
  group_by(Sample) %>%
  summarise(total_rel_abundance = sum(relative_abundance, na.rm = TRUE))
print("Sums of relative abundances at species level after filtering")
print(temp)



meta_families <- meta_filtered %>%
  # Extract only rows that contain a family-level clade
  filter(str_detect(Species, "f__[^|]+")) %>%
  # Extract the family name and keep relevant columns
  transmute(
    Sample,
    Family = str_extract(Species, "f__[^|]+"),
    relative_abundance,
    est_read_counts = estimated_number_of_reads_from_the_clade
  ) %>%
  # Aggregate by sample and family
  group_by(Sample, Family) %>%
  summarise(
    relative_abundance = sum(relative_abundance, na.rm = TRUE),
    est_read_counts = sum(est_read_counts, na.rm = TRUE),
    .groups = "drop"
  )

# confirm that family-level rel abundances sum to 100
temp <- meta_families %>%
  group_by(Sample) %>%
  summarise(total_rel_abundance = sum(relative_abundance, na.rm = TRUE))
print("Sums of family-level relative abundances")
print(temp)

# -----
# Family-level bar plot
# -----
# Compute mean abundance per family
family_means <- meta_families %>%
  complete(Sample, Family, fill = list(relative_abundance = 0)) %>%  # missing families = 0
  group_by(Family) %>%
  summarise(mean_rel_abund = mean(relative_abundance), .groups = "drop")

# Get the top 15 families
top_families <- family_means %>%
  arrange(desc(mean_rel_abund)) %>%
  slice_head(n = 15) %>%
  pull(Family)


# Collapse all non-top families into "Other"
meta_families_plot <- meta_families %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%
  group_by(Sample, Family) %>%
  summarise(
    relative_abundance = sum(relative_abundance),
    est_read_counts = sum(est_read_counts),
    .groups = "drop"
  )

# Recompute mean abundance after collapsing (including "Other")
family_means <- meta_families_plot %>%
  group_by(Family) %>%
  summarise(mean_rel_abund = mean(relative_abundance, na.rm = TRUE)) %>%
  ungroup()

# Add family labels with average %
meta_families_plot <- meta_families_plot %>%
  left_join(family_means, by = "Family") %>%
  mutate(family_label = paste0(Family, " (", round(mean_rel_abund, 1), "%)"))

# Order legend by mean abundance
# Top families by mean abundance, Other at the end
meta_families_plot <- meta_families_plot %>%
  mutate(family_label = fct_reorder(family_label, mean_rel_abund, .desc = TRUE))

# Convert family_label to factor if it's not already
meta_families_plot$family_label <- factor(meta_families_plot$family_label)

# Detect the level that starts with "Other"
other_level <- levels(meta_families_plot$family_label)[grepl("^Other", levels(meta_families_plot$family_label))]

# Relevel Other to be last
meta_families_plot$family_label <- fct_relevel(meta_families_plot$family_label, other_level, after = Inf)

meta_families_plot <- meta_families_plot %>%
  group_by(Sample) %>%
  mutate(relative_abundance = relative_abundance / sum(relative_abundance) * 100) %>%
  ungroup()

colors <- c(
  "#0074D9",
  "#FF851B",
  "#2ECC40",
  "#FF4136",
  "#B10DC9",
  "#FFDC00",
  "#39CCCC",
  "#85144B",
  "#01FF70",
  "#001F3F",
  "#F012BE",
  "#3D9970",
  "#7FDBFF",
  "#A52A2A",
  "#FF7F50",
  "#808080"
  )

# Plot stacked barplot
p <- ggplot(meta_families_plot, aes(x = Sample, y = relative_abundance, fill = family_label)) +
  geom_bar(stat = "identity") +
  # scale_fill_viridis_d(option = "turbo") +  # discrete turbo palette
  scale_fill_manual(values = colors) +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Family (avg %)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "families_summary_stacked_barplot.pdf", plot = p, width = 20, height = 6)


# reorder the samples by family similarity
meta_wide <- meta_families_plot %>%
  select(Sample, family_label, relative_abundance) %>%
  pivot_wider(names_from = family_label, values_from = relative_abundance, values_fill = 0)

# Remove Sample column
abund_matrix <- as.matrix(meta_wide[,-1])
rownames(abund_matrix) <- meta_wide$Sample

# Compute Bray-Curtis distance
dist_matrix <- vegdist(abund_matrix, method = "bray")

# Hierarchical clustering
hc <- hclust(dist_matrix, method = "average")

# Extract sample order
# Extract ordered sample names
sample_order <- rownames(abund_matrix)[hc$order]

# Reorder factor in plotting data
meta_families_plot <- meta_families_plot %>%
  mutate(Sample = factor(Sample, levels = sample_order))

p2 <- ggplot(meta_families_plot, aes(x = Sample, y = relative_abundance, fill = family_label)) +
  geom_bar(stat = "identity") +
  # scale_fill_viridis_d(option = "turbo") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "families_summary_stacked_barplot_clustered.pdf", plot = p2, width = 20, height = 6)

# -----




# Align with metadata
meta_filtered$Sample <- gsub("_profile", "", meta_filtered$Sample)
sample_order <- unique(meta_filtered$Sample)
meta_df <- meta_df %>%
  filter(!is.na(File_ID) & File_ID %in% sample_order) %>%
  arrange(match(File_ID, sample_order))
print(all(meta_df$File_ID %in% sample_order))

# Filter metadata for important covariates
meta_df <- meta_df[, c(
  "record_id",
  "Block",
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
  "File_ID"
)]
meta_df$smke <- ifelse(meta_df$smke %in% c("1", "2", "3"), 1,
                       ifelse(meta_df$smke == "0", 0, NA))
meta_df$alco <- ifelse(meta_df$alco %in% c("1", "2"), 0,
                       ifelse(meta_df$alco %in% c("3", "4"), 1, NA))

# Run alpha and beta diversity
otu_mat <- meta_filtered %>%
  select(Sample, Species, relative_abundance) %>%
  pivot_wider(
    names_from = Sample,
    values_from = relative_abundance,
    values_fill = 0
  )

taxa <- otu_mat$Species
otu_rel <- as.matrix(otu_mat[, -1])
rownames(otu_rel) <- taxa

otu <- otu_table(otu_rel, taxa_are_rows = TRUE)

# just keep taxonomy as one string for now
tax <- tax_table(matrix(taxa, ncol = 1))
rownames(tax) <- taxa
colnames(tax) <- "FullTaxonomy"

ps <- phyloseq(otu, tax)
# match metadata
all(colnames(otu_rel) %in% rownames(meta_df))
# Create sample_data object
sdata <- sample_data(meta_df)
# Add to existing phyloseq object
ps <- merge_phyloseq(ps, sdata)

alpha_div <- estimate_richness(ps, measures = "Shannon")
head(alpha_div)

meta_df$alpha_diversity <- alpha_div[rownames(meta_df), "Shannon"]

# add total reads and number of species
temp <- meta_filtered %>%
  group_by(Sample) %>%
  summarise(total_reads = sum(estimated_number_of_reads_from_the_clade, na.rm = TRUE),
            num_species = n_distinct(Species)) %>%
  as.data.frame()
rownames(temp) <- temp$Sample
meta_df$est_total_reads <- temp[rownames(meta_df), "total_reads"]
meta_df$num_species <- temp[rownames(meta_df), "num_species"]

# Does alpha diversity depend on input read count
meta_df$Block <- factor(meta_df$Block)
lm_fit <- lm(alpha_diversity ~ est_total_reads, data = meta_df)
intercept <- coef(lm_fit)[1]
slope <- coef(lm_fit)[2]
r2 <- summary(lm_fit)$r.squared
eq_label <- paste0("y = ", round(slope, 4), "x + ", round(intercept, 2))
r2_label <- paste0("R² = ", round(r2, 3))
p <- ggplot(meta_df, aes(x = est_total_reads, y = alpha_diversity)) +
  geom_point(aes(color = Block), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +   # linear regression line
  annotate("text", x = max(meta_df$est_total_reads)*0.6, 
           y = max(meta_df$alpha_diversity)*0.9, 
           label = eq_label, color = "black", size = 4) +
  annotate("text", x = max(meta_df$est_total_reads)*0.6, 
           y = max(meta_df$alpha_diversity)*0.85, 
           label = r2_label, color = "black", size = 4) +
  theme_minimal()
ggsave("alpha_div_vs_total_reads_block.pdf", plot = p, width = 6, height = 6, dpi = 300)

p <- ggplot(meta_df, aes(x = est_total_reads, y = alpha_diversity)) +
  geom_point(aes(color = sarc_status), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +   # linear regression line
  annotate("text", x = max(meta_df$est_total_reads)*0.6, 
          y = max(meta_df$alpha_diversity)*0.9, 
          label = eq_label, color = "black", size = 4) +
  annotate("text", x = max(meta_df$est_total_reads)*0.6, 
           y = max(meta_df$alpha_diversity)*0.85, 
           label = r2_label, color = "black", size = 4) +
  theme_minimal()
ggsave("alpha_div_vs_total_reads_sarc.pdf", plot = p, width = 6, height = 6, dpi = 300)

lm_fit <- lm(num_species ~ est_total_reads, data = meta_df)
intercept <- coef(lm_fit)[1]
slope <- coef(lm_fit)[2]
r2 <- summary(lm_fit)$r.squared
eq_label <- paste0("y = ", round(slope, 4), "x + ", round(intercept, 2))
r2_label <- paste0("R² = ", round(r2, 3))
p <- ggplot(meta_df, aes(x = est_total_reads, y = num_species)) +
  geom_point(aes(color = Block), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +   # linear regression line
  annotate("text", x = max(meta_df$est_total_reads)*0.6, 
          y = max(meta_df$num_species)*0.9, 
          label = eq_label, color = "black", size = 4) +
  annotate("text", x = max(meta_df$est_total_reads)*0.6, 
           y = max(meta_df$num_species)*0.85, 
           label = r2_label, color = "black", size = 4) +
  theme_minimal()
ggsave("num_species_vs_total_reads_block.pdf", plot = p, width = 6, height = 6, dpi = 300)

p <- ggplot(meta_df, aes(x = est_total_reads, y = num_species, color = sarc_status)) +
  geom_point(aes(color = sarc_status), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +   # linear regression line
  annotate("text", x = max(meta_df$est_total_reads)*0.6, 
          y = max(meta_df$num_species)*0.9, 
          label = eq_label, color = "black", size = 4) +
  annotate("text", x = max(meta_df$est_total_reads)*0.6, 
           y = max(meta_df$num_species)*0.85, 
           label = r2_label, color = "black", size = 4) +
  theme_minimal()
ggsave("num_species_vs_total_reads_sarc.pdf", plot = p, width = 6, height = 6, dpi = 300)


# beta diversity
meta_df$sarc_status <- as.numeric(meta_df$sarc_status)
meta_df$age_def <- as.numeric(meta_df$age_def)
ps <- phyloseq(otu, tax)
# match metadata
all(colnames(otu_rel) %in% rownames(meta_df))
# Create sample_data object
sdata <- sample_data(meta_df)
# Add to existing phyloseq object
ps <- merge_phyloseq(ps, sdata)

bray_dist <- phyloseq::distance(ps, method = "bray")

ord <- ordinate(ps, method = "PCoA", distance = bray_dist)

p <- plot_ordination(ps, ord, color = "sarc_status") + 
  geom_point(size = 3) +
  theme_minimal()
ggsave("beta_diversity_sarc_status.pdf", plot = p, width = 6, height = 6, dpi = 300)
p <- plot_ordination(ps, ord, color = "Block") + 
  geom_point(size = 3) +
  theme_minimal()
ggsave("beta_diversity_block.pdf", plot = p, width = 6, height = 6, dpi = 300)
p <- plot_ordination(ps, ord, color = "est_total_reads") + 
  geom_point(size = 3) +
  theme_minimal()
ggsave("beta_diversity_total_reads.pdf", plot = p, width = 6, height = 6, dpi = 300)


meta <- data.frame(sample_data(ps))
# keep only samples with non-NA Treatment
meta_clean <- meta[!is.na(meta$sarc_status), ]
# subset the distance matrix to the same samples
bray_clean <- as.matrix(bray_dist)[rownames(meta_clean), rownames(meta_clean)]
bray_clean <- as.dist(bray_clean)
adonis_res <- adonis2(bray_clean ~ sarc_status, data = meta_clean, permutations = 999)
adonis_res
# adonis2(formula = bray_clean ~ sarc_status, data = meta_clean, permutations = 999)
#           Df SumOfSqs      R2      F Pr(>F)   
# Model      1    0.531 0.01429 2.8423  0.005 **
# Residual 196   36.582 0.98571                  
# Total    197   37.112 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

meta <- data.frame(sample_data(ps))
# keep only samples with non-NA Treatment
meta_clean <- meta[!is.na(meta$Block), ]
# subset the distance matrix to the same samples
bray_clean <- as.matrix(bray_dist)[rownames(meta_clean), rownames(meta_clean)]
bray_clean <- as.dist(bray_clean)
adonis_res <- adonis2(bray_clean ~ Block, data = meta_clean, permutations = 999)
adonis_res
# adonis2(formula = bray_clean ~ Block, data = meta_clean, permutations = 999)
#           Df SumOfSqs      R2     F Pr(>F)
# Model      6    1.302 0.03438 1.157  0.194
# Residual 195   36.579 0.96562             
# Total    201   37.881 1.00000   

meta <- data.frame(sample_data(ps))
# keep only samples with non-NA Treatment
meta_clean <- meta[!is.na(meta$est_total_reads), ]
# subset the distance matrix to the same samples
bray_clean <- as.matrix(bray_dist)[rownames(meta_clean), rownames(meta_clean)]
bray_clean <- as.dist(bray_clean)
adonis_res <- adonis2(bray_clean ~ est_total_reads, data = meta_clean, permutations = 999)
adonis_res
# adonis2(formula = bray_clean ~ est_total_reads, data = meta_clean, permutations = 999)
#            Df SumOfSqs      R2      F Pr(>F)    
# Model      1    2.595 0.06851 14.711  0.001 ***
# Residual 200   35.286 0.93149                  
# Total    201   37.881 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




# DESeq2
meta_df$age_def_scaled <- (meta_df$age_def - mean(meta_df$age_def, na.rm = TRUE)) / sd(meta_df$age_def, na.rm = TRUE)
meta_df$nutr_score <- as.numeric(meta_df$nutr_score)
meta_df$nutr_score_scaled <- (meta_df$nutr_score - mean(meta_df$nutr_score, na.rm = TRUE)) / sd(meta_df$nutr_score, na.rm = TRUE)
meta_df$bmi <- as.numeric(meta_df$bmi)
meta_df$bmi_scaled <- (meta_df$bmi - mean(meta_df$bmi, na.rm = TRUE)) / sd(meta_df$bmi, na.rm = TRUE)
meta_df$sex <- ifelse(meta_df$sex %in% c("0","1"), meta_df$sex, NA)
meta_df$sex <- factor(meta_df$sex, levels = c("0","1"))

cor(meta_df$nutr_score_scaled, meta_df$bmi_scaled, method = "spearman", use = "complete.obs")
# 0.065

# remake the OTU table with estimated read counts
otu_mat <- meta_filtered %>%
  select(Sample, Species, estimated_number_of_reads_from_the_clade) %>%
  pivot_wider(
    names_from = Sample,
    values_from = estimated_number_of_reads_from_the_clade,
    values_fill = 0
  )

taxa <- otu_mat$Species
otu_rel <- as.matrix(otu_mat[, -1])
rownames(otu_rel) <- taxa

otu <- otu_table(otu_rel, taxa_are_rows = TRUE)

# just keep taxonomy as one string for now
tax <- tax_table(matrix(taxa, ncol = 1))
rownames(tax) <- taxa
colnames(tax) <- "FullTaxonomy"

ps <- phyloseq(otu, tax)
# match metadata
all(colnames(otu_rel) %in% rownames(meta_df))
# Create sample_data object
sdata <- sample_data(meta_df)
# Add to existing phyloseq object
ps <- merge_phyloseq(ps, sdata)

# No NAs in DESeq2
ps_clean <- subset_samples(ps, 
               !is.na(age_def_scaled) &
               !is.na(sex) &
               !is.na(smke) &
               !is.na(alco) &
               !is.na(nutr_score_scaled) &
               !is.na(bmi_scaled) &
               !is.na(sarc_status))
dds <- phyloseq_to_deseq2(ps_clean, ~ age_def_scaled + sex + smke + alco + nutr_score_scaled + bmi_scaled + sarc_status)
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res)
write.csv(res_df, "deseq2_results.csv", row.names = TRUE)

# Visualize DESeq2 results in dotplot
# Keep only significant species
sig_res <- res_df %>%
  filter(!is.na(padj), padj < 0.05)
sig_res <- sig_res %>%
  mutate(
    genus_species = str_extract(rownames(sig_res), "g__[^|]+\\|s__[^|]+"),  # extract genus|species
    genus_species = str_replace_all(genus_species, "g__|s__", "")          # remove prefixes
  )

# Create the plot
p <- ggplot(sig_res, aes(
    x = log2FoldChange,
    y = reorder(genus_species, log2FoldChange),
    color = -log10(pvalue),
    size = baseMean
  )) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_viridis_c(option = "turbo") +
  labs(
    x = "Log2 Fold Change",
    y = "Species",
    color = "-log10(p-value)",
    size = "baseMean",
    title = "Differential abundance of significant species"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))
ggsave("deseq2_dotplot.pdf", plot = p, width = 10, height = 10)




# Write meta_filtered
write.csv(meta_filtered, "meta_filtered.csv", row.names = FALSE)

# Write meta_df
write.csv(meta_df, "meta_df.csv", row.names = FALSE)

