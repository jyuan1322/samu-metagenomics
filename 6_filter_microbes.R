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

input_dir <- "/data/local/jy1008/SaMu/metaphlan_out_09172025/all_merged_fastqs"

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
meta_df <- do.call(rbind, lapply(files, read_metaphlan))

# Convert relative abundance to numeric
meta_df <- meta_df %>%
  mutate(relative_abundance = as.numeric(relative_abundance))

# Keep only species-level groups
meta_species <- meta_df %>%
  filter(str_detect(clade_name, "s__")) %>%
  mutate(
    # Extract species name
    species = str_extract(clade_name, "s__[^|]+")
  ) %>%
  group_by(Sample, species) %>%
  summarise(
    relative_abundance = sum(relative_abundance, na.rm = TRUE),
    estimated_number_of_reads_from_the_clade = sum(estimated_number_of_reads_from_the_clade, na.rm = TRUE),
    .groups = "drop"
  )


# Elbow plot
# Define a sequence of abundance thresholds, e.g., 0% to 1%
abundance_thresholds <- seq(0, 0.1, by = 0.001)  # 0% to 10%
prevalence_thresholds <- c(1, 5, 10, 20)  # number of samples a species must appear in

# Build elbow data for multiple prevalence thresholds
elbow_df <- expand_grid(
  Threshold = abundance_thresholds,
  MinSamples = prevalence_thresholds
) %>%
  mutate(
    CladesRetained = map2_int(Threshold, MinSamples, ~ {
      meta_species %>%
        group_by(species) %>%
        summarize(n_samples_above = sum(relative_abundance >= .x), .groups = "drop") %>%
        filter(n_samples_above >= .y) %>%
        nrow()
    })
  )

# Plot multiple lines
p <- ggplot(elbow_df, aes(x = Threshold * 100, y = CladesRetained, color = factor(MinSamples))) +
  geom_line(size = 1.2) +
  theme_minimal() +
  labs(
    x = "Relative Abundance Threshold (%)",
    y = "Number of Species Passing Filter",
    color = "Min Samples Passing Threshold",
    title = "Elbow Plot for Abundance Filtering"
  )
ggsave("elbow_plot_species_level_09232025.pdf", plot = p, width = 10, height = 6)


# ----------



# Identify species that pass the threshold in at least 5 samples
species_keep <- meta_species %>%
  group_by(species) %>%
  summarize(n_samples_above_threshold = sum(relative_abundance >= 0.1)) %>%
  filter(n_samples_above_threshold >= 5) %>%
  pull(species)

# Filter original dataframe
meta_filtered <- meta_df %>%
  filter(clade_name %in% clade_keep)

# Print number of clades before and after filtering
n_before <- meta_df %>%
  pull(clade_name) %>%
  unique() %>%
    length()
n_after <- meta_filtered %>%
  pull(clade_name) %>%
  unique() %>%
    length()
cat("Number of clades before filtering:", n_before, "\n")
cat("Number of clades after filtering:", n_after, "\n")
# Number of clades before filtering: 2128 
# Number of clades after filtering: 785


# meta_wide <- meta_filtered %>%
#   select(Sample, clade_name, relative_abundance) %>%
#   pivot_wider(names_from = Sample, values_from = relative_abundance, values_fill = 0)

meta_families <- meta_filtered %>%
  filter(str_detect(clade_name, "f__[^|]+$")) %>%   # only entries ending with f__something
  transmute(
    Sample,
    family = str_extract(clade_name, "f__[^|]+$"),
    relative_abundance,
    est_read_counts = estimated_number_of_reads_from_the_clade
  )

# -----
# Family-level bar plot
# -----
# Compute mean abundance per family
family_means <- meta_families %>%
  complete(Sample, family, fill = list(relative_abundance = 0)) %>%  # missing families = 0
  group_by(family) %>%
  summarise(mean_rel_abund = mean(relative_abundance), .groups = "drop")

# Get the top 15 families
top_families <- family_means %>%
  arrange(desc(mean_rel_abund)) %>%
  slice_head(n = 15) %>%
  pull(family)





# Collapse all non-top families into "Other"
meta_families_plot <- meta_families %>%
  mutate(family = ifelse(family %in% top_families, family, "Other")) %>%
  group_by(Sample, family) %>%
  summarise(
    relative_abundance = sum(relative_abundance),
    est_read_counts = sum(est_read_counts),
    .groups = "drop"
  )

# Recompute mean abundance after collapsing (including "Other")
family_means <- meta_families_plot %>%
  group_by(family) %>%
  summarise(mean_rel_abund = mean(relative_abundance, na.rm = TRUE)) %>%
  ungroup()

# Add family labels with average %
meta_families_plot <- meta_families_plot %>%
  left_join(family_means, by = "family") %>%
  mutate(family_label = paste0(family, " (", round(mean_rel_abund, 1), "%)"))

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

# Plot stacked barplot
p <- ggplot(meta_families_plot, aes(x = Sample, y = relative_abundance, fill = family_label)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "turbo") +  # discrete turbo palette
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
  scale_fill_viridis_d(option = "turbo") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "families_summary_stacked_barplot_clustered.pdf", plot = p2, width = 20, height = 6)

