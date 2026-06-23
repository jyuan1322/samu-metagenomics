# =============================================================================
# 03_differential_abundance.R
#
# DESeq2 differential abundance on species-level estimated read counts,
# adjusting for covariates, with apeglm shrinkage and a faceted dotplot of
# significant species. Consumes the intermediates from 01_load_and_filter.R
# (and the diversity-augmented meta_df from 02_diversity.R if present).
#
# Run:  Rscript 03_differential_abundance.R
# =============================================================================
source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  library(phyloseq); library(DESeq2); library(tidytext)
})

setwd(OUTPUT_DIR)

meta_filtered <- readRDS(META_FILTERED_RDS)
meta_df       <- readRDS(META_DF_RDS)

meta_filtered$Sample <- gsub("_profile", "", meta_filtered$Sample)

# ---------------------------------------------------------------------------
# Scale continuous covariates and set factor levels
# ---------------------------------------------------------------------------
for (cv in SCALE_COVARIATES) {
  meta_df[[cv]] <- as.numeric(meta_df[[cv]])
  meta_df[[paste0(cv, "_scaled")]] <-
    (meta_df[[cv]] - mean(meta_df[[cv]], na.rm = TRUE)) /
    sd(meta_df[[cv]], na.rm = TRUE)
}
meta_df$sex <- factor(meta_df$sex, levels = c("M", "F"))

# ---------------------------------------------------------------------------
# Build a species-level count matrix (estimated reads) aligned to metadata
# ---------------------------------------------------------------------------
otu_mat <- meta_filtered %>%
  select(Sample, Species, estimated_number_of_reads_from_the_clade) %>%
  pivot_wider(names_from = Sample,
              values_from = estimated_number_of_reads_from_the_clade,
              values_fill = 0)
taxa <- otu_mat$Species
otu_counts <- as.matrix(otu_mat[, -1]); rownames(otu_counts) <- taxa

rownames(meta_df) <- meta_df$File_ID
common <- intersect(colnames(otu_counts), rownames(meta_df))
message(sprintf("Samples shared between counts and metadata: %d", length(common)))
otu_counts <- otu_counts[, common]
meta_df    <- meta_df[common, ]
stopifnot(all(colnames(otu_counts) == rownames(meta_df)))

otu <- otu_table(otu_counts, taxa_are_rows = TRUE)
tax <- tax_table(matrix(taxa, ncol = 1))
rownames(tax) <- taxa; colnames(tax) <- "FullTaxonomy"
ps <- merge_phyloseq(phyloseq(otu, tax), sample_data(meta_df))

# ---------------------------------------------------------------------------
# DESeq2 (drop any sample with NA in a model term)
# ---------------------------------------------------------------------------
model_terms <- all.vars(DESEQ_DESIGN)
complete_samples <- rownames(meta_df)[complete.cases(meta_df[, model_terms])]
ps_clean <- prune_samples(complete_samples, ps)
message(sprintf("Samples entering DESeq2: %d", length(complete_samples)))

dds <- phyloseq_to_deseq2(ps_clean, DESEQ_DESIGN)
dds <- DESeq(dds)

res <- lfcShrink(dds, coef = DESEQ_COEF, type = "apeglm")
res_df <- as.data.frame(res)
write.csv(res_df, tag_filename("deseq2_results.csv"), row.names = TRUE)
saveRDS(dds, tag_filename("deseq2_dds.rds"))

# ---------------------------------------------------------------------------
# Mean relative abundance in the reference group (for annotation/joins)
# ---------------------------------------------------------------------------
ps_rel  <- transform_sample_counts(ps, function(x) x / sum(x))
ps_ref  <- subset_samples(ps_rel, get(GROUP_VAR) == GROUP_LEVELS[1])
otu_ref <- otu_table(ps_ref)
otu_ref <- if (taxa_are_rows(ps_ref)) as.matrix(otu_ref) else t(as.matrix(otu_ref))
mean_rel_ref <- data.frame(Species = rownames(otu_ref),
                           mean_rel_ref = rowMeans(otu_ref))

# ---------------------------------------------------------------------------
# Dotplot of significant species, faceted by |log2FC|
# ---------------------------------------------------------------------------
sig_res <- res_df %>%
  filter(!is.na(padj), padj < PADJ_CUTOFF) %>%
  mutate(genus_species = str_extract(rownames(.), "s__[^|]+"),
         genus_species = str_replace_all(genus_species, "g__|s__", ""),
         genus_species = str_replace_all(genus_species, "_", " "),
         fullname = rownames(.)) %>%
  left_join(mean_rel_ref, by = c("fullname" = "Species"))

if (nrow(sig_res) > 0) {
  sig_res$highFC <- ifelse(abs(sig_res$log2FoldChange) >= LFC_FACET_CUTOFF,
                           paste0("log2 Fold Change >= ", LFC_FACET_CUTOFF),
                           paste0("log2 Fold Change < ", LFC_FACET_CUTOFF))
  sig_res$genus_species <- reorder_within(sig_res$genus_species,
                                          sig_res$log2FoldChange,
                                          sig_res$highFC)

  p <- ggplot(sig_res, aes(log2FoldChange, genus_species,
                           color = -log10(pvalue))) +
    geom_point(size = 4) +
    geom_vline(xintercept = 0, color = "black") +
    scale_color_viridis_c(option = "plasma") +
    facet_wrap(~highFC, nrow = 1, scales = "free_y") +
    scale_y_reordered() +
    labs(x = "Log2 Fold Change", y = "Species",
         color = "-log10(p-value)",
         title = "Differential abundance of significant species") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  ggsave(tag_filename("deseq2_dotplot_split.pdf"), p, width = 12, height = 6)
  ggsave(tag_filename("deseq2_dotplot_split.png"), p,
         width = 12, height = 6, dpi = 300)
} else {
  message("No species passed padj < ", PADJ_CUTOFF, "; skipping dotplot.")
}

write.csv(meta_filtered, tag_filename("meta_filtered.csv"), row.names = FALSE)
write.csv(meta_df,       tag_filename("meta_df_FullSaMu.csv"), row.names = FALSE)
message("03_differential_abundance.R complete.")
