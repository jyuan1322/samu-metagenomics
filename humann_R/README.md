# HUMAnN Downstream Analysis (R)

Sibling of `metagenomics_R/`, but for HUMAnN pathway/gene family output
(from `7_merge_humann.sh`) instead of MetaPhlAn taxonomic profiles. No
alpha/beta diversity — not really meaningful for a functional profile.
Filtering adds a coefficient-of-variation step on top of abundance/
prevalence; the endpoint is MaAsLin3 instead of DESeq2.

**This folder does not rebuild File_ID or recode covariates.** It reads the
already-recoded, complete-case-filtered `meta_df.rds` written by
`metagenomics_R/01_load_and_filter.R` (via `METAGENOMICS_META_DF_RDS` in
`config.R`), so that logic — `build_file_id()`, `recode_metadata()`,
`FILE_ID_BLOCKS` — lives in exactly one place. Run `metagenomics_R`'s
pipeline (at least through `01_load_and_filter.R`) for the cohort/settings
you want before running this one.

## Scripts

| File | Role |
|------|------|
| `config.R` | Paths (including `METAGENOMICS_META_DF_RDS`), filtering thresholds, sample-ID matching regexes, MaAsLin3 model spec. |
| `utils.R` | `read_humann_table`, `extract_sample_id`, `diagnose_sample_matching`, `split_unstratified`, `feature_stats`, `filter_by_abundance_cv`, `tag_filename` (same convention as `metagenomics_R/utils.R`). |
| `01_load_and_filter.R` | Load the joined+renormalized HUMAnN table, keep unstratified rows, align sample columns to `meta_df$File_ID`, filter by abundance + prevalence + CV, save diagnostic plot + filtered intermediates. |
| `02_maaslin.R` | Run MaAsLin3 with `sarc_status_bin` (already factored NoSarc/Sarc from `meta_df.rds`) as the fixed effect. |

Run in order from this directory, after `metagenomics_R/01_load_and_filter.R`
has produced `meta_df.rds`:

```r
Rscript 01_load_and_filter.R
Rscript 02_maaslin.R
```

## Sample-ID matching — read this before your first run

`humann_join_tables` names each merged-table column after the input
filename minus its `_<file_name>.tsv` suffix — e.g. a file named
`SaMu9_host_removed_R1R2_combined_humann_pathabundance.tsv` becomes a column
like `SaMu9_host_removed_R1R2_combined_humann` (possibly `_Abundance`
appended). That is **not** the same string as `File_ID` (`"SaMu9"`) — both
MetaPhlAn and HUMAnN run on the same fastq per sample, so the same File_ID is
embedded as a *prefix* either way, but MetaPhlAn's output filename convention
is exactly `<File_ID>_profile.txt` (File_ID = the whole prefix), while
HUMAnN's carries extra pipeline-stage suffixes (`_host_removed_R1R2_combined_humann`,
from `3_remove_host_reads.sh` / `1_combine_L001_L002.sh` / `6_humann.sh`), so
File_ID is only a leading substring there.

`extract_sample_id()` resolves this in two steps:

1. `STRIP_SUFFIX_REGEX` removes known trailing unit suffixes (`_Abundance`, `_RPK`, `_CPM`).
2. `FILE_ID_REGEXES` tries each block's extraction pattern in turn against
   what's left — the *same* patterns as `FILE_ID_BLOCKS`' regex values in
   `metagenomics_R/config.R`, just applied to the HUMAnN column name instead
   of a metadata column value. This reuses the same block coverage
   `build_file_id()` already handles on the MetaPhlAn side, rather than
   assuming every sample is block1/SaMu-prefixed. Order matters where one
   pattern is a prefix of another (block2's `\d+_\d+` is a prefix of
   block5/89a's patterns, so it's listed last) — keep `FILE_ID_REGEXES` in
   sync with `FILE_ID_BLOCKS` if that list changes.

`01_load_and_filter.R` calls `diagnose_sample_matching()` and prints match
counts plus a few unmatched examples from each side before failing loudly —
check that output on your first run rather than trusting a nonzero result
count blindly.

## Why a CV filter, and why on top of abundance

HUMAnN pathway tables usually contain a handful of near-universal,
high-abundance pathways (core carbon/energy metabolism, etc.) present at
similar levels in essentially every sample. These pass a standard abundance/
prevalence filter easily but carry little information for a between-group
association test — there's nothing to associate, since everyone has about
the same amount. `MIN_CV` drops these directly: features must vary (sd/mean
above threshold, computed over nonzero values only) in addition to being
abundant and prevalent enough to trust.

`01_load_and_filter.R` writes `abundance_vs_cv.png` before applying the
filter — the diagnostic to look at when tuning `MIN_MEAN_ABUNDANCE` and
`MIN_CV`, same role as the elbow plot in `metagenomics_R/01_load_and_filter.R`.

## MaAsLin3, not MaAsLin2

Current lab-supported successor to MaAsLin2; tests both abundance and
prevalence associations and handles compositionality better in general use.

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biobakery/maaslin3")
```

If you add `RANDOM_EFFECTS` (e.g. `"record_id"` for repeated sampling),
random effects behave differently in MaAsLin3 than MaAsLin2 — check the
`maaslin3` manual's random-effects section and sanity-check fixed-effect
results with and without the random effect before trusting them.

## Configuration notes

- **`METAGENOMICS_META_DF_RDS`** — must point at the `meta_df.rds` from the
  specific `metagenomics_R` run (i.e. `OUTPUT_DIR`, `RUN_TAG`,
  `EXTREME_CASES_ONLY`) you want this analysis aligned to.
- **`FEATURE_TABLE`** — `"pathabundance"` or `"genefamilies"`. Start with
  pathways; gene families are far larger and slower.
- **`MIN_MEAN_ABUNDANCE` / `MIN_PREVALENCE` / `MIN_CV`** — start permissive,
  tighten using `abundance_vs_cv.png` and the resulting feature count.
- **`MAASLIN_NORMALIZATION`** — defaults to `"TSS"` (MaAsLin3's own
  recommended default). Since the input is already relab-normalized and then
  feature-filtered, TSS re-normalizes relative to the filtered feature set,
  not the original whole-sample total — standard practice, but worth
  remembering when comparing effect sizes across filtering runs. Set to
  `"NONE"` to keep the original relab scale instead.

## Known gaps / things to verify on first run

- `feature_stats()` computes mean/CV over **nonzero** values only, not all
  samples — `MIN_MEAN_ABUNDANCE` is "mean when present," not "mean across
  the whole cohort." `MIN_PREVALENCE` is what excludes low-prevalence,
  high-when-present features.
- `SAMPLE_ID_EXTRACT_REGEX` is confirmed to match the SaMu-prefixed naming
  seen so far, but hasn't been checked against every sequencing block — see
  the matching section above.


# e.g. 
# python plot_maaslin_results.py
#   --input /data/local/jy1008/SaMu/results/latest/humann_R/maaslin3_pathabundance_07162026/all_results.tsv \
#   --metadata sarc_status_bin \
#   --omit-file omit_pathways.txt \
#   --label-top-n 3 \
#   --output maaslin_volcano_top_features


python plot_maaslin_results.py \
    --input /data/local/jy1008/SaMu/results/latest/humann_R/maaslin3_pathabundance_07162026/all_results.tsv \
    --metadata sarc_status_bin \
    --top-n 15 --label-top-n 3 \
    --omit-file humann_omit_pathways.txt \
    --abundance-table /data/local/jy1008/SaMu/results/latest/humann_R/filtered_features_abundance_07162026.csv \
    --sample-metadata /data/local/jy1008/SaMu/results/latest/humann_R/sample_metadata_for_heatmap_07162026.csv \
    --sample-id-col File_ID --group-col sarc_status_bin \
    --output maaslin_volcano_top_features