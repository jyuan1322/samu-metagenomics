# =============================================================================
# 02_diversity.R
#
# Family- and genus-level stacked barplots (Bray-Curtis clustered), alpha
# diversity (normalized Shannon) vs read depth and by group, and beta
# diversity (PCoA + PERMANOVA). Consumes the intermediates from
# 01_load_and_filter.R.
#
# Run:  Rscript 02_diversity.R
# =============================================================================
source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
  library(forcats); library(vegan); library(phyloseq)
  library(tidytext); library(rstatix)
})

setwd(OUTPUT_DIR)

meta_filtered <- readRDS(META_FILTERED_RDS)
meta_df       <- readRDS(META_DF_RDS)

# ---------------------------------------------------------------------------
# Helper: collapse a taxonomic rank to top taxa + "Other", then a
# Bray-Curtis-clustered stacked barplot. Used for both family and genus.
# ---------------------------------------------------------------------------
rank_barplot <- function(df, rank_regex, rank_name, out_stem) {
  # Pull the rank label out of each species string and aggregate.
  ranked <- df %>%
    filter(str_detect(Species, rank_regex)) %>%
    transmute(Sample,
              Rank = str_extract(Species, rank_regex),
              relative_abundance,
              est_read_counts = estimated_number_of_reads_from_the_clade) %>%
    group_by(Sample, Rank) %>%
    summarise(relative_abundance = sum(relative_abundance, na.rm = TRUE),
              est_read_counts    = sum(est_read_counts, na.rm = TRUE),
              .groups = "drop")

  write.csv(ranked, tag_filename(paste0(out_stem, "_abundances.csv")),
            row.names = FALSE)

  # Top taxa: mean abundance over the cutoff, excluding "unclassified".
  rank_means <- ranked %>%
    complete(Sample, Rank, fill = list(relative_abundance = 0)) %>%
    group_by(Rank) %>%
    summarise(mean_rel_abund = mean(relative_abundance), .groups = "drop")
  top_taxa <- rank_means %>%
    filter(mean_rel_abund > BARPLOT_MEAN_ABUND_CUTOFF,
           !str_detect(Rank, "unclassified")) %>%
    pull(Rank)

  plot_df <- ranked %>%
    mutate(Rank = ifelse(Rank %in% top_taxa, Rank, "Other")) %>%
    group_by(Sample, Rank) %>%
    summarise(relative_abundance = sum(relative_abundance),
              est_read_counts    = sum(est_read_counts), .groups = "drop")

  # Label with mean %, order legend by abundance, push "Other" last.
  plot_means <- plot_df %>%
    group_by(Rank) %>%
    summarise(mean_rel_abund = mean(relative_abundance, na.rm = TRUE))
  plot_df <- plot_df %>%
    left_join(plot_means, by = "Rank") %>%
    mutate(rank_label = paste0(Rank, " (", round(mean_rel_abund, 1), "%)"),
           rank_label = fct_reorder(rank_label, mean_rel_abund, .desc = TRUE))
  other_lvl <- levels(plot_df$rank_label)[grepl("^Other", levels(plot_df$rank_label))]
  plot_df$rank_label <- fct_relevel(plot_df$rank_label, other_lvl, after = Inf)

  plot_df <- renormalize_relab(plot_df)

  # Order samples by Bray-Curtis similarity (hierarchical clustering).
  wide <- plot_df %>%
    select(Sample, rank_label, relative_abundance) %>%
    pivot_wider(names_from = rank_label, values_from = relative_abundance,
                values_fill = 0)
  m <- as.matrix(wide[, -1]); rownames(m) <- wide$Sample
  hc <- hclust(vegdist(m, method = "bray"), method = "average")
  plot_df$Sample <- factor(plot_df$Sample, levels = rownames(m)[hc$order])

  p <- ggplot(plot_df, aes(Sample, relative_abundance, fill = rank_label)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = BARPLOT_COLORS,
                      name = paste0(rank_name, " (avg %)")) +
    labs(x = "Sample", y = "Relative Abundance (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(tag_filename(paste0(out_stem, "_stacked_barplot_clustered.pdf")),
         p, width = 20, height = 6)
  ggsave(tag_filename(paste0(out_stem, "_stacked_barplot_clustered.png")),
         p, width = 20, height = 6, dpi = 300)
}

rank_barplot(meta_filtered, "f__[^|]+", "Family", "families_summary")
rank_barplot(meta_filtered, "g__[^|]+", "Genus",  "genera_summary")

# ---------------------------------------------------------------------------
# Build a species-level phyloseq object (relative abundances) for diversity
# ---------------------------------------------------------------------------
meta_filtered$Sample <- gsub("_profile", "", meta_filtered$Sample)
sample_order <- unique(meta_filtered$Sample)
meta_df <- meta_df %>%
  filter(!is.na(File_ID) & File_ID %in% sample_order) %>%
  arrange(match(File_ID, sample_order))
stopifnot(all(meta_df$File_ID %in% sample_order))

otu_mat <- meta_filtered %>%
  select(Sample, Species, relative_abundance) %>%
  pivot_wider(names_from = Sample, values_from = relative_abundance,
              values_fill = 0)
taxa <- otu_mat$Species
otu_rel <- as.matrix(otu_mat[, -1]); rownames(otu_rel) <- taxa

otu <- otu_table(otu_rel, taxa_are_rows = TRUE)
tax <- tax_table(matrix(taxa, ncol = 1));
rownames(tax) <- taxa; colnames(tax) <- "FullTaxonomy"

rownames(meta_df) <- meta_df$File_ID
ps <- merge_phyloseq(phyloseq(otu, tax), sample_data(meta_df))

# ---------------------------------------------------------------------------
# Alpha diversity (normalized Shannon)
# ---------------------------------------------------------------------------
alpha_div <- estimate_richness(ps, measures = c("Shannon"))
alpha_div$Shannon_norm <- alpha_div$Shannon / log(nrow(tax_table(ps)))

meta_df$alpha_diversity      <- alpha_div[rownames(meta_df), "Shannon"]
meta_df$alpha_diversity_norm <- alpha_div[rownames(meta_df), "Shannon_norm"]

read_stats <- meta_filtered %>%
  group_by(Sample) %>%
  summarise(total_reads = sum(estimated_number_of_reads_from_the_clade,
                              na.rm = TRUE),
            num_species = n_distinct(Species)) %>%
  as.data.frame()
rownames(read_stats) <- read_stats$Sample
meta_df$est_total_reads <- read_stats[rownames(meta_df), "total_reads"]
meta_df$num_species     <- read_stats[rownames(meta_df), "num_species"]

message("Sum of estimated reads from clades: ", sum(meta_df$est_total_reads))
message("Median estimated reads from clades: ", median(meta_df$est_total_reads))

# Alpha diversity vs read depth.
lm_fit <- lm(alpha_diversity_norm ~ est_total_reads, data = meta_df)
r2_label <- paste0("R^2 = ", round(summary(lm_fit)$r.squared, 3))
p <- ggplot(meta_df, aes(est_total_reads, alpha_diversity_norm)) +
  geom_point(aes(color = .data[[GROUP_VAR]]), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = max(meta_df$est_total_reads) * 0.6,
           y = max(meta_df$alpha_diversity_norm) * 0.85,
           label = r2_label, color = "black", size = 4) +
  theme_minimal() +
  labs(x = "Summed Total Reads from Clade",
       y = "Normalized Shannon Alpha Diversity")
ggsave(tag_filename("alpha_div_normed_vs_total_reads.pdf"), p,
       width = 6, height = 6)
ggsave(tag_filename("alpha_div_normed_vs_total_reads.png"), p,
       width = 6, height = 6, dpi = 300)

# Alpha diversity by group.
p <- ggplot(meta_df, aes(.data[[GROUP_VAR]], alpha_diversity_norm)) +
  geom_boxplot(aes(fill = .data[[GROUP_VAR]]),
               outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(x = "Group", y = "Normalized Shannon Alpha Diversity")
ggsave(tag_filename("alpha_diversity_normed_by_group.pdf"), p,
       width = 6, height = 6)
ggsave(tag_filename("alpha_diversity_normed_by_group.png"), p,
       width = 6, height = 6, dpi = 300)

write.csv(meta_df, tag_filename("metadata_with_alpha_diversity.csv"),
          row.names = FALSE)

# Wilcoxon test of alpha diversity between groups (BH-adjusted).
group_formula <- as.formula(paste("alpha_diversity_norm ~", GROUP_VAR))
alpha_wilcox <- meta_df %>%
  wilcox_test(group_formula, ref.group = GROUP_LEVELS[1]) %>%
  adjust_pvalue(method = "BH")
write.csv(alpha_wilcox,
          tag_filename("alpha_diversity_normed_group_wilcox.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Beta diversity (Bray-Curtis PCoA + PERMANOVA)
# ---------------------------------------------------------------------------
bray_dist <- phyloseq::distance(ps, method = "bray")
ord <- ordinate(ps, method = "PCoA", distance = bray_dist)

p <- plot_ordination(ps, ord, color = GROUP_VAR) +
  geom_point(size = 3) + theme_minimal()
ggsave(tag_filename("beta_diversity_group.pdf"), p, width = 6, height = 6)
ggsave(tag_filename("beta_diversity_group.png"), p,
       width = 6, height = 6, dpi = 300)

meta <- data.frame(sample_data(ps))
meta_clean <- meta[!is.na(meta[[GROUP_VAR]]), ]
bray_clean <- as.dist(as.matrix(bray_dist)[rownames(meta_clean),
                                           rownames(meta_clean)])
adonis_formula <- as.formula(paste("bray_clean ~", GROUP_VAR))
adonis_res <- adonis2(adonis_formula, data = meta_clean, permutations = 9999)
print(adonis_res)
write.csv(as.data.frame(adonis_res),
          tag_filename("adonis_group_results.csv"), row.names = TRUE)

saveRDS(meta_df, META_DF_RDS)   # updated with diversity columns
message("02_diversity.R complete.")
