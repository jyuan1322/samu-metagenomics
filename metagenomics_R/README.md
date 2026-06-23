# Metagenomics Downstream Analysis (R)

Takes MetaPhlAn per-sample profiles (the `*_profile.txt` output of the
preprocessing pipeline), joins them to sample metadata, filters, and produces
diversity and differential-abundance results. All study-specific values live
in `config.R`; the numbered scripts should not need editing between studies.

## Scripts

| File | Role |
|------|------|
| `config.R` | All paths, thresholds, the File_ID block mapping, and the DESeq2 design. Edit this per study. |
| `utils.R` | Shared helpers: `build_file_id`, `recode_metadata`, `read_metaphlan`, `tag_filename`, `renormalize_relab`, palette. |
| `01_load_and_filter.R` | Metadata recode + File_ID link, load profiles, species-level collapse, abundance/prevalence filter, elbow plot. Saves intermediates. |
| `02_diversity.R` | Family/genus stacked barplots (Bray-Curtis clustered), alpha diversity (normalized Shannon) + Wilcoxon, beta diversity (PCoA + PERMANOVA). |
| `03_differential_abundance.R` | DESeq2 on species-level counts with covariate adjustment + apeglm shrinkage, faceted dotplot of significant species. |

Run them in order from this directory:

```r
Rscript 01_load_and_filter.R
Rscript 02_diversity.R
Rscript 03_differential_abundance.R
```

State passes between scripts via two `.rds` files written into `OUTPUT_DIR`
(`meta_filtered.rds`, `meta_df.rds`), so each script can be re-run on its own
once `01` has been run.

## Configuration notes

- **`INPUT_DIR` / `OUTPUT_DIR` / `METADATA_CSV`** — point these at the study.
  `OUTPUT_DIR` receives every plot, CSV, and `.rds`.
- **`RUN_TAG`** — a suffix (e.g. a date) appended to versioned outputs via
  `tag_filename()`, so a new run does not overwrite a previous one. Set to `""`
  to disable.
- **`FILE_ID_BLOCKS`** — the per-block column-to-regex mapping that builds the
  `File_ID` linking metadata rows to profile files. Each sequencing block
  stores its file reference in a different column with a different ID format.
  Blocks are applied in order; a non-empty match in a later block overrides an
  earlier one for the same row. Add/remove entries as blocks change — this
  replaces the previously hardcoded four-block sequence.
- **Filtering** — `ABUND_THRESHOLD` (%) and `MIN_SAMPLES` set the
  abundance/prevalence filter; the elbow plot from `01` is the diagnostic for
  choosing them. `BARPLOT_MEAN_ABUND_CUTOFF` controls the top-taxa cutoff in
  the stacked barplots.
- **Model** — `DESEQ_DESIGN` is the design formula (put the grouping variable
  last), `DESEQ_COEF` is the coefficient handed to `lfcShrink` and must match a
  `resultsNames(dds)` entry, and `SCALE_COVARIATES` lists the continuous
  covariates that get z-scored into `<name>_scaled` columns referenced by the
  formula. `GROUP_VAR` / `GROUP_LEVELS` define the binary contrast used across
  all three scripts.

## Behavior changes from the original single script

- Family and genus barplots are produced by one `rank_barplot()` helper rather
  than two near-duplicate blocks; output filenames are unchanged in spirit
  (`families_summary_*`, `genera_summary_*`, with `RUN_TAG` appended).
- The four hardcoded File_ID blocks are now data in `config.R`.
- Output filenames use `tag_filename()` instead of hardcoded dates.
- The DESeq2 complete-case filter now derives its columns from the design
  formula (`all.vars(DESEQ_DESIGN)`) instead of a hand-listed `subset_samples`
  call, so changing the formula automatically changes which samples are kept.
- Group-specific names (`sarc_status_bin`, `NoSarc`) are replaced by
  `GROUP_VAR` / `GROUP_LEVELS` throughout; plots/CSVs say "group" generically.

## Requirements

R packages: `dplyr`, `tidyr`, `purrr`, `readr`, `ggplot2`, `stringr`,
`forcats`, `tools`, `viridis`, `vegan`, `phyloseq`, `DESeq2`, `tidytext`,
`rstatix`, and `apeglm` (for `lfcShrink`). `phyloseq`, `DESeq2`, and `apeglm`
are from Bioconductor.
