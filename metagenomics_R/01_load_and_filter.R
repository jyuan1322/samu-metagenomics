# =============================================================================
# 01_load_and_filter.R
#
# Load metadata and MetaPhlAn profiles, build the File_ID link, recode
# covariates, filter to the analysis cohort, collapse to species level, apply
# the abundance/prevalence filter, renormalize, and save the intermediates
# consumed by 02_diversity.R and 03_differential_abundance.R.
#
# Run from this directory:  Rscript 01_load_and_filter.R
# =============================================================================
source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(readr); library(ggplot2); library(stringr)
})

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
setwd(OUTPUT_DIR)

# ---------------------------------------------------------------------------
# Metadata: load, build File_ID, recode, cohort filter
# ---------------------------------------------------------------------------
meta_df <- read.csv(METADATA_CSV, header = TRUE, stringsAsFactors = FALSE)

meta_df <- build_file_id(meta_df)            # uses FILE_ID_BLOCKS from config
rownames(meta_df) <- meta_df[[RECORD_ID_COL]]

# Restrict to the cohort flag (e.g. Full.SaMu == 1), if configured.
if (!is.na(FULLSAMU_COL)) {
  meta_df <- meta_df[meta_df[[FULLSAMU_COL]] == 1, ]
}

meta_df <- recode_metadata(meta_df, extreme_cases_only = EXTREME_CASES_ONLY)          # smke/alco/sex/sarc_status recode

# check where NAs are produced
sapply(meta_df[c("sarc_status","age_def","bmi","nutr_score")],
       function(x) sum(is.na(suppressWarnings(as.numeric(as.character(x))))))

write.csv(meta_df,
          tag_filename("metadata_FullSaMu_recoded.csv"),
          row.names = FALSE)

# Drop subjects with no sequencing File_ID (present for other assays only).
meta_df <- meta_df[!is.na(meta_df$File_ID), ]

# Keep the covariates of interest, then require them to be complete.
meta_df <- meta_df[, META_KEEP_COLS]
meta_df <- meta_df[complete.cases(meta_df[, COLS_OF_INTEREST]), ]
message(sprintf("Samples after cohort + complete-case filtering: %d",
                nrow(meta_df)))
print("total subjects in each class:")
table(meta_df$sarc_status)

# ---------------------------------------------------------------------------
# MetaPhlAn profiles: load only the files for retained subjects
# ---------------------------------------------------------------------------
files <- list.files(INPUT_DIR, pattern = "_profile.txt$", full.names = TRUE)
files <- files[basename(files) %in% paste0(meta_df$File_ID, "_profile.txt")]

comb_df <- do.call(rbind, lapply(files, read_metaphlan)) %>%
  mutate(relative_abundance = as.numeric(relative_abundance))

# Species-level entries only (drop strain "|t__" rows).
meta_species <- comb_df %>%
  filter(str_detect(clade_name, "s__")) %>%
  filter(!str_detect(clade_name, "\\|t__")) %>%
  rename(Species = clade_name) %>%
  group_by(Sample, Species) %>%
  summarise(
    relative_abundance = sum(relative_abundance, na.rm = TRUE),
    estimated_number_of_reads_from_the_clade =
      sum(estimated_number_of_reads_from_the_clade, na.rm = TRUE),
    .groups = "drop"
  )

message("Per-sample species-level relative abundance sums (expect ~100):")
print(meta_species %>% group_by(Sample) %>%
        summarise(total = sum(relative_abundance, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# Elbow plot: clades retained vs abundance threshold, by prevalence cutoff
# ---------------------------------------------------------------------------
elbow_df <- expand_grid(
  threshold   = ELBOW_ABUND_THRESHOLDS,
  min_samples = ELBOW_PREVALENCE_THRESHOLDS
) %>%
  mutate(CladesRetained = map2_int(threshold, min_samples, ~ {
    meta_species %>%
      group_by(Species) %>%
      summarise(n_samples_above = sum(relative_abundance >= .x),
                .groups = "drop") %>%
      filter(n_samples_above >= .y) %>%
      nrow()
  }))

p_elbow <- ggplot(elbow_df,
                  aes(threshold, CladesRetained,
                      color = factor(min_samples))) +
  geom_line(size = 1.2) +
  scale_x_continuous(breaks = seq(0, 2, by = 0.1)) +
  theme_minimal() +
  labs(x = "Relative Abundance Threshold (%)",
       y = "Number of Species Passing Filter",
       color = "Min Samples Passing Threshold",
       title = "Elbow Plot for Abundance Filtering")
ggsave(tag_filename("elbow_plot_species_level.pdf"), p_elbow,
       width = 10, height = 6)
ggsave(tag_filename("elbow_plot_species_level.png"), p_elbow,
       width = 10, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Apply the abundance/prevalence filter, then renormalize to 100%
# ---------------------------------------------------------------------------
species_keep <- meta_species %>%
  group_by(Species) %>%
  summarise(n_above = sum(relative_abundance >= ABUND_THRESHOLD)) %>%
  filter(n_above >= MIN_SAMPLES) %>%
  pull(Species)

meta_filtered <- meta_species %>% filter(Species %in% species_keep)

cat("Number of clades before filtering:",
    length(unique(meta_species$Species)), "\n")
cat("Number of clades after filtering:",
    length(unique(meta_filtered$Species)), "\n")

meta_filtered <- renormalize_relab(meta_filtered)

message("Per-sample sums after filtering + renormalization (expect ~100):")
print(meta_filtered %>% group_by(Sample) %>%
        summarise(total = sum(relative_abundance, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# Save intermediates for the downstream scripts
# ---------------------------------------------------------------------------
write.csv(meta_filtered, tag_filename("meta_filtered.csv"), row.names = FALSE)
saveRDS(meta_filtered, META_FILTERED_RDS)
saveRDS(meta_df,       META_DF_RDS)

message("01_load_and_filter.R complete.")
