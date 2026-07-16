# =============================================================================
# 03_pathway_hierarchy.R
#
# Parse the BioCyc SmartTable export of chained "Ontology - direct parents of
# entity" columns (built manually on biocyc.org, scoped to MetaCyc, from
# surviving_pathway_ids.txt) and print each surviving pathway's ancestor
# hierarchy, truncated at the domain-meaningful "Pathways" root (dropping the
# generic Pathway-Tools framework classes above it: Generalized-Reactions,
# FRAMES, THINGS).
#
# IMPORTANT LIMITATION: once a pathway has more than one direct parent (a
# "//"-separated cell), BioCyc's chained export loses the exact
# parent<->grandparent pairing — column N+1 is the union of direct-parents-
# of-everything-in-column-N for that row, not a disambiguated per-branch
# lineage. This script therefore treats each column as a "breadth layer"
# (everything N steps up from the pathway), not a fully resolved tree. That's
# sufficient for a hierarchy printout, but not for anything requiring exact
# edges (e.g. a GOBPPARENTS-style igraph distance calculation, as used for
# GO in run_ORA.R).
#
# Input: a tab-delimited export of the SmartTable you built on biocyc.org —
# column 1 = pathway ID, columns 2+ = repeated "Ontology - direct parents of
# entity" transform columns, each cell possibly containing multiple values
# separated by " // ". Export via the SmartTable's download/export option
# and set PATHWAY_ONTOLOGY_EXPORT_TSV in config.R to point at the file.
#
# Outputs (in OUTPUT_DIR):
#   pathway_hierarchy.txt              — human-readable indented hierarchy per pathway
#   pathway_hierarchy.csv              — tidy long format: pathway_id, level, ancestor
#   redundancy_structural_candidates.csv — pathways with "Super-Pathways" at level 1
#                                           (composite pathways; check against SUB_PATHWAYS_EXPORT_TSV
#                                           if set, otherwise these are candidates to check manually)
#   redundancy_statistical_pairs.csv   — pathway pairs with |correlation| >= REDUNDANCY_COR_THRESHOLD
#                                           on the filtered abundance matrix
#   redundancy_cor_histogram.png/pdf   — distribution of pairwise correlations, to sanity-check
#                                           REDUNDANCY_COR_THRESHOLD isn't arbitrary for this data
#
# Run from this directory, after 01_load_and_filter.R:  Rscript 03_pathway_hierarchy.R
# =============================================================================
source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(purrr); library(readr); library(ggplot2)
})

setwd(OUTPUT_DIR)

# ---------------------------------------------------------------------------
# Load the BioCyc export
# ---------------------------------------------------------------------------
lineage_raw <- read.delim(PATHWAY_ONTOLOGY_EXPORT_TSV, header = TRUE,
                           check.names = FALSE, stringsAsFactors = FALSE,
                           sep = "\t")
names(lineage_raw)[1] <- "pathway_id"
level_cols <- names(lineage_raw)[-1]
names(lineage_raw)[-1] <- paste0("level", seq_along(level_cols))
level_cols <- names(lineage_raw)[-1]

message(sprintf("Loaded ancestor export: %d pathways x %d transform columns",
                nrow(lineage_raw), length(level_cols)))

# ---------------------------------------------------------------------------
# Cross-check against the surviving pathway list — flag anything missing so
# a partial/stale export doesn't silently produce an incomplete report.
# ---------------------------------------------------------------------------
surviving_ids <- readLines(tag_filename("surviving_pathway_ids_only.txt"))
missing_from_export <- setdiff(surviving_ids, lineage_raw$pathway_id)
extra_in_export <- setdiff(lineage_raw$pathway_id, surviving_ids)
if (length(missing_from_export) > 0) {
  message(sprintf("WARNING: %d surviving pathway(s) missing from the BioCyc export: %s",
                  length(missing_from_export), paste(head(missing_from_export, 10), collapse = ", ")))
}
if (length(extra_in_export) > 0) {
  message(sprintf("NOTE: %d pathway(s) in the export are not in surviving_pathway_ids_only.txt (ignored): %s",
                  length(extra_in_export), paste(head(extra_in_export, 10), collapse = ", ")))
}

# ---------------------------------------------------------------------------
# Parse each row into ordered "breadth layers" of ancestor classes, splitting
# multi-parent cells on " // " and truncating once ROOT_CLASS_NAME (default
# "Pathways") is reached in a layer — everything past that point is the
# generic Pathway-Tools frame hierarchy (Generalized-Reactions/FRAMES/THINGS),
# not a meaningful pathway category.
# ---------------------------------------------------------------------------
parse_lineage_row <- function(row_values, root_class = ROOT_CLASS_NAME) {
  layers <- list()
  for (v in row_values) {
    if (is.na(v) || trimws(v) == "") break
    items <- str_split(v, "\\s*//\\s*")[[1]]
    items <- unique(items[items != ""])
    layers[[length(layers) + 1]] <- items
    if (root_class %in% items) break
  }
  layers
}

pathway_lineages <- setNames(
  lapply(seq_len(nrow(lineage_raw)), function(i) {
    parse_lineage_row(as.character(lineage_raw[i, level_cols]))
  }),
  lineage_raw$pathway_id
)

n_unresolved <- sum(!map_lgl(pathway_lineages, ~ ROOT_CLASS_NAME %in% unlist(.x)))
if (n_unresolved > 0) {
  message(sprintf(
    "WARNING: %d pathway(s) never reached '%s' within %d transform columns — chain further on biocyc.org and re-export.",
    n_unresolved, ROOT_CLASS_NAME, length(level_cols)))
}

# ---------------------------------------------------------------------------
# Human-readable indented report
# ---------------------------------------------------------------------------
hierarchy_lines <- character()
for (pwy in names(pathway_lineages)) {
  hierarchy_lines <- c(hierarchy_lines, pwy)
  layers <- pathway_lineages[[pwy]]
  for (i in seq_along(layers)) {
    hierarchy_lines <- c(hierarchy_lines,
                         sprintf("%s-> %s", strrep("  ", i), paste(layers[[i]], collapse = " | ")))
  }
  hierarchy_lines <- c(hierarchy_lines, "")
}
writeLines(hierarchy_lines, tag_filename("pathway_hierarchy.txt"))

# ---------------------------------------------------------------------------
# Tidy long-format CSV: one row per (pathway_id, level, ancestor class) —
# useful for programmatic filtering/grouping later (e.g. "which pathways
# share a level-1 ancestor of Cofactor-Biosynthesis").
# ---------------------------------------------------------------------------
tidy_rows <- imap_dfr(pathway_lineages, function(layers, pwy) {
  if (length(layers) == 0) return(tibble(pathway_id = character(), level = integer(), ancestor = character()))
  imap_dfr(layers, function(items, lvl) {
    tibble(pathway_id = pwy, level = lvl, ancestor = items)
  })
})
write_csv(tidy_rows, tag_filename("pathway_hierarchy.csv"))

# ---------------------------------------------------------------------------
# Redundancy check 1 (structural): pathways whose level-1 ancestor list
# includes "Super-Pathways" are themselves composite — their reaction list
# includes one or more sub-pathways' reactions wholesale, so if both a
# superpathway and its own sub-pathway survive filtering, they aren't
# independent evidence.
#
# This only flags CANDIDATES from data already in hand (the hierarchy just
# parsed). Confirming which surviving pathway(s) each superpathway actually
# overlaps with requires a separate BioCyc lookup: add a "Sub-Pathways"
# transform column to a SmartTable of just these candidate IDs, export it,
# and set SUB_PATHWAYS_EXPORT_TSV in config.R — if set, this section
# cross-checks it automatically; if not, the candidate list alone is written
# out for manual follow-up.
# ---------------------------------------------------------------------------
superpathway_candidates <- tidy_rows %>%
  filter(level == 1, ancestor == "Super-Pathways") %>%
  distinct(pathway_id) %>%
  pull(pathway_id)

message(sprintf("%d surviving pathway(s) flagged as superpathways (structural redundancy candidates)",
                length(superpathway_candidates)))
write_csv(tibble(pathway_id = superpathway_candidates),
          tag_filename("redundancy_structural_candidates.csv"))

if (!is.null(SUB_PATHWAYS_EXPORT_TSV) && !is.na(SUB_PATHWAYS_EXPORT_TSV) && file.exists(SUB_PATHWAYS_EXPORT_TSV)) {
  sub_raw <- read.delim(SUB_PATHWAYS_EXPORT_TSV, header = TRUE, check.names = FALSE,
                        stringsAsFactors = FALSE, sep = "\t")
  names(sub_raw)[1] <- "pathway_id"
  names(sub_raw)[2] <- "sub_pathways"

  structural_overlaps <- sub_raw %>%
    mutate(sub_pathways = str_split(sub_pathways, "\\s*//\\s*")) %>%
    tidyr::unnest(sub_pathways) %>%
    filter(sub_pathways %in% surviving_ids) %>%
    rename(superpathway = pathway_id, sub_pathway = sub_pathways)

  message(sprintf("%d confirmed superpathway/sub-pathway overlap(s) among surviving pathways",
                  nrow(structural_overlaps)))
  write_csv(structural_overlaps, tag_filename("redundancy_structural_confirmed.csv"))
} else {
  message("SUB_PATHWAYS_EXPORT_TSV not set/found — skipping structural cross-check; ",
          "see redundancy_structural_candidates.csv for pathways to check manually.")
}

# ---------------------------------------------------------------------------
# Redundancy check 2 (statistical): pathway pairs that behave nearly
# identically across your samples, regardless of formal MetaCyc structure —
# e.g. both driven by the same dominant organism. Spearman rather than
# Pearson since relative-abundance data is compositional/skewed. Computed on
# the FILTERED feature matrix (post 01_load_and_filter.R), not the raw table.
# ---------------------------------------------------------------------------
filtered_mat <- readRDS(FILTERED_FEATURES_RDS)
cor_mat <- cor(t(filtered_mat), method = "spearman", use = "pairwise.complete.obs")

cor_vals <- cor_mat[upper.tri(cor_mat)]
p_cor_hist <- ggplot(tibble(cor = cor_vals), aes(cor)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  geom_vline(xintercept = c(-1, 1) * REDUNDANCY_COR_THRESHOLD, color = "firebrick", linetype = "dashed") +
  theme_minimal() +
  labs(x = "Pairwise Spearman correlation (filtered pathways)",
       y = "Number of pathway pairs",
       title = "Distribution of pairwise pathway correlations")
ggsave(tag_filename("redundancy_cor_histogram.pdf"), p_cor_hist, width = 8, height = 6)
ggsave(tag_filename("redundancy_cor_histogram.png"), p_cor_hist, width = 8, height = 6, dpi = 300)

high_cor_idx <- which(abs(cor_mat) >= REDUNDANCY_COR_THRESHOLD & upper.tri(cor_mat), arr.ind = TRUE)
redundant_pairs <- tibble(
  pathway_1 = rownames(cor_mat)[high_cor_idx[, 1]],
  pathway_2 = colnames(cor_mat)[high_cor_idx[, 2]],
  correlation = cor_mat[high_cor_idx]
) %>%
  arrange(desc(abs(correlation)))

message(sprintf("%d pathway pair(s) with |Spearman correlation| >= %.2f",
                nrow(redundant_pairs), REDUNDANCY_COR_THRESHOLD))
write_csv(redundant_pairs, tag_filename("redundancy_statistical_pairs.csv"))

message(sprintf("03_pathway_hierarchy.R complete. %d pathways parsed. Output in: %s",
                length(pathway_lineages), OUTPUT_DIR))