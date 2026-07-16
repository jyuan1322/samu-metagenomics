# =============================================================================
# config.R — central configuration for the HUMAnN downstream analysis
# (joined pathway/gene family tables -> metadata join -> filtering -> MaAsLin3).
#
# Source this at the top of every numbered script:
#   source("config.R")
#
# Sibling of metagenomics_R/config.R, same conventions (RUN_TAG, tag_filename,
# GROUP_VAR/GROUP_LEVELS). Deliberately does NOT re-derive File_ID or recode
# covariates — it reuses the already-recoded, complete-case-filtered
# meta_df.rds written by metagenomics_R/01_load_and_filter.R, so the File_ID
# construction and smke/alco/sex/sarc_status_bin recoding logic lives in
# exactly one place.
# =============================================================================

# ---------------------------------------------------------------------------
# Paths — inputs
# ---------------------------------------------------------------------------
# Directory produced by 7_merge_humann.sh (contains joined_*_relab.tsv etc.)
HUMANN_MERGED_DIR <- "/data/local/jy1008/SaMu/humann_out_07152026/all_merged_fastqs_merged"

PATHABUNDANCE_FILE <- file.path(HUMANN_MERGED_DIR, "joined_pathabundance_relab.tsv")
GENEFAMILIES_FILE  <- file.path(HUMANN_MERGED_DIR, "joined_genefamilies_relab.tsv")

# The metagenomics_R run whose meta_df.rds should be reused here. Must match
# the OUTPUT_DIR of whichever metagenomics_R/01_load_and_filter.R run you
# want this analysis aligned to (same cohort filter, same EXTREME_CASES_ONLY
# setting, same recoded sarc_status_bin factor levels).
METAGENOMICS_OUTPUT_DIR <- "/data/local/jy1008/SaMu/results/latest/metagenomics_R_test2"
METAGENOMICS_META_DF_RDS <- file.path(METAGENOMICS_OUTPUT_DIR, "meta_df.rds")

# ---------------------------------------------------------------------------
# Paths — outputs
# ---------------------------------------------------------------------------
OUTPUT_DIR <- "/data/local/jy1008/SaMu/results/latest/humann_R"

# Same convention as metagenomics_R: appended before the extension via
# tag_filename() so a new run does not overwrite a previous one.
# RUN_TAG <- "07162026"
RUN_TAG <- "07162026"

# Intermediates passed between the numbered scripts (written into OUTPUT_DIR,
# untagged — matches how META_DF_RDS is saved untagged in metagenomics_R).
FILTERED_FEATURES_RDS <- "filtered_features.rds"
META_ALIGNED_RDS       <- "meta_aligned.rds"

# ---------------------------------------------------------------------------
# Which table to analyze
# ---------------------------------------------------------------------------
# "pathabundance" or "genefamilies". Pathways are the natural first pass;
# gene families are ~2 orders of magnitude larger (see 7_merge_humann.sh
# output file sizes) and slower to load/filter — only switch once the
# pathway-level result is in hand.
FEATURE_TABLE <- "pathabundance"

# ---------------------------------------------------------------------------
# Sample-ID matching: HUMAnN column name -> File_ID
# ---------------------------------------------------------------------------
# humann_join_tables names each sample column after the input filename minus
# its "_<file_name>.tsv" suffix (e.g. a file named
# "SaMu9_host_removed_R1R2_combined_humann_pathabundance.tsv" becomes column
# "SaMu9_host_removed_R1R2_combined_humann" [+ "_Abundance" for some HUMAnN
# versions/units]) — this is NOT the same string as File_ID
# ("SaMu9"), which is why simple suffix-stripping isn't enough on its own.
#
# Both tools run on the same fastq per sample, so the same File_ID is embedded
# as a prefix in both MetaPhlAn's and HUMAnN's output names — but the exact
# *string* differs: MetaPhlAn's output is exactly "<File_ID>_profile.txt", so
# File_ID is the whole prefix, nothing else. HUMAnN's input filename carries
# extra pipeline-stage suffixes from upstream steps (e.g.
# "SaMu9_host_removed_R1R2_combined_humann_pathabundance.tsv" — the
# "_host_removed_R1R2_combined_humann" part comes from 3_remove_host_reads.sh /
# 1_combine_L001_L002.sh / 6_humann.sh naming), so File_ID is only a prefix
# substring, not the whole thing.
#
# STRIP_SUFFIX_REGEX removes known trailing unit suffixes first (_Abundance
# etc.); FILE_ID_REGEXES then tries each block's extraction pattern in turn
# against what's left, same patterns as FILE_ID_BLOCKS in
# metagenomics_R/config.R (just applied to a different source string — the
# HUMAnN column name rather than a metadata column value). This reuses the
# same block coverage build_file_id() already handles on the MetaPhlAn side,
# rather than assuming every sample is block1/SaMu-prefixed.
#
# Keep these regexes in sync with FILE_ID_BLOCKS' regex values if that list
# changes. Verify against real column names before trusting alignment counts:
#   source("config.R"); source("utils.R")
#   colnames(read_humann_table(PATHABUNDANCE_FILE))
# and compare against meta_df$File_ID from METAGENOMICS_META_DF_RDS.
STRIP_SUFFIX_REGEX <- "_Abundance$|_RPK$|_CPM$"
# Order matters: more specific patterns must come first, since a less
# specific pattern that's a *substring* of a more specific one will match
# first and silently truncate the extraction otherwise. Concretely: block1's
# "SaMu[0-9]+" matches "SaMu165" as a substring buried inside a block5-style
# name like "4915_01_121_SaMu165_host_removed..." — if tried before block5's
# pattern, it grabs just "SaMu165" instead of the full "4915_01_121_SaMu165",
# which then fails to match meta_df$File_ID. So: block89a and block5 (both
# contain "SaMu\d+" as a trailing substring of a longer pattern) must be
# tried before block1's bare "SaMu[0-9]+"; block2's "\d+_\d+" is a prefix of
# both numeric patterns and must be tried last of all.
FILE_ID_REGEXES <- c(
  "\\d+_\\d+_Libr\\d+_SaMu\\d+",    # block89a (most specific)
  "\\d+_\\d+_\\d+_SaMu\\d+",        # block5
  "SaMu[0-9]+",                     # block1 (substring of block5/89a — must follow them)
  "\\d+_\\d+"                       # block2 (prefix of block5/89a numeric parts — least specific, tried last)
)

# ---------------------------------------------------------------------------
# Filtering thresholds
# ---------------------------------------------------------------------------
# Minimum mean relative abundance (fraction, not %; HUMAnN relab sums to 1.0
# per sample), computed over samples where the feature is nonzero.
MIN_MEAN_ABUNDANCE <- 1e-4

# Minimum prevalence: fraction of samples in which the feature is nonzero.
MIN_PREVALENCE <- 0.10

# Minimum coefficient of variation (sd/mean over nonzero values). Drops
# pathways/gene families that are abundant in nearly every sample but don't
# vary between samples — high abundance, low information for an association
# test. Tune against abundance_vs_cv.png from 01_load_and_filter.R.
MIN_CV <- 0.1

# ---------------------------------------------------------------------------
# Elbow plot sweeps (diagnostic for choosing the three thresholds above) —
# same idea as metagenomics_R's ELBOW_ABUND_THRESHOLDS / ELBOW_PREVALENCE_THRESHOLDS:
# sweep a threshold on the x-axis, one line per minimum-prevalence requirement,
# showing how many features survive. Two sweeps here (abundance, CV) instead
# of metagenomics_R's one, since this pipeline filters on both.
# ---------------------------------------------------------------------------
# Mean relative abundance thresholds to sweep. Log-spaced (not linear like
# metagenomics_R's seq(0, 2, by=0.02)) since HUMAnN relab spans several
# orders of magnitude rather than metagenomics_R's ~0-2% range.
ELBOW_ABUND_THRESHOLDS <- 10 ^ seq(-6, -1, length.out = 60)

# Coefficient-of-variation thresholds to sweep.
ELBOW_CV_THRESHOLDS <- seq(0, 3, by = 0.05)

# Minimum-prevalence *fractions* compared as separate lines on both elbow
# plots (converted internally to a sample count: ceiling(fraction * n_samples)).
# Fractions rather than metagenomics_R's raw counts (c(1,5,10,20)) since this
# config's own MIN_PREVALENCE is already fraction-based.
ELBOW_PREVALENCE_FRACTIONS <- c(0, 0.05, 0.10, 0.20, 0.50)

# ---------------------------------------------------------------------------
# MaAsLin3 model
# ---------------------------------------------------------------------------
# Matches GROUP_VAR in metagenomics_R/config.R; meta_df$sarc_status_bin
# arrives already recoded and factored (levels NoSarc, Sarc) from
# METAGENOMICS_META_DF_RDS, so no re-coercion should be necessary — 02_maaslin.R
# checks this and only coerces as a fallback.
FIXED_EFFECTS  <- c("sarc_status_bin")
RANDOM_EFFECTS <- c()   # e.g. c("record_id") if you have repeated measures

# TSS+LOG are MaAsLin3's own recommended/validated defaults. Since the input
# is already HUMAnN relab-normalized and then feature-filtered, TSS
# re-normalizes relative to the filtered feature set, not the original
# whole-sample total — standard practice, but worth remembering when
# interpreting effect sizes. Set to "NONE" to keep the original relab scale.
MAASLIN_NORMALIZATION    <- "TSS"
MAASLIN_TRANSFORM        <- "LOG"
MAASLIN_MAX_SIGNIFICANCE <- 0.1