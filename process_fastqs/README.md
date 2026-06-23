# Metagenomics Preprocessing Pipeline

Shotgun metagenomics preprocessing, from raw sequencer tar archives to
taxonomic profiles. The pipeline runs in numbered stages; each stage has a
plain-loop version (runs serially in your shell) and, where the original
workflow had one, a separate SLURM version. All paths, database locations,
resource requests, and environment names live in `config.sh` — edit that one
file to point the pipeline at a new study.

## Stages

| Stage | Plain-loop script | SLURM script | What it does |
|-------|-------------------|--------------|--------------|
| 0 | `0_collect_tar_fastqs.sh` | `0_submit_extract_jobs.sh` (+ `0_extract_fastq.sh` worker) | Extract `fastq.gz` from raw tar archives |
| 1 | `1_combine_L001_L002.sh` | — | Concatenate L001/L002 lanes per sample |
| 1 (multi) | `1_combine_L001_L002_multiple_replicates.sh` | — | Same, pooling across replicate directories |
| 2 | `2_qc.sh` | `2_qc_array.slurm` (+ `make_qc_list.sh`) | fastp trimming + fastqc + multiqc |
| 3 | `3_remove_host_reads.sh` | `3_remove_host_reads_job_submit.sh` (+ `3_remove_host_reads_job.sh` worker) | bowtie2 host alignment, keep unmapped pairs |
| 4 | `4_kraken.sh` | — | Kraken2 + Bracken classification |
| 5 | `5_metaphlan.sh` | `5_submit_metaphlan_jobs.sh` | MetaPhlAn profiling |

`config.sh` and `common.sh` are sourced by the scripts and are not run directly.

## Setup

1. Edit `config.sh`:
   - `BASE_DIR` — project data root. All working directories (`raw/`,
     `mapping/`, `host_removed/`, `kraken_out/`, `metaphlan_out/`, `logs/`)
     are derived from it.
   - `TAR_FILES` — list of raw tar archives (stage 0 only).
   - `SUBFOLD` — optional per-batch subdirectory used by stages 3 and 5.
     Set to `""` to work directly on the top-level directories.
   - Reference DBs: `HOST_BOWTIE2_INDEX`, `KRAKEN_DB`, `METAPHLAN_BOWTIE2_DB`.
   - `CPUS`, the `SLURM_*` resource values, and the micromamba env names
     (`METAGEN_ENV`, `METAPHLAN_ENV`).

2. For local runs, activate the relevant environment yourself first
   (e.g. `micromamba activate metagen-env`). SLURM scripts activate the
   environment themselves via `common.sh`'s `activate_env`.

## Running

Local (serial) example, top to bottom:

```bash
# Stage 0 — extract archives listed in config.sh
./0_collect_tar_fastqs.sh

# Stage 1 — combine lanes (single dir)
./1_combine_L001_L002.sh "$BASE_DIR/raw/<run>" "$BASE_DIR/raw/all_merged_fastqs"
#   ...or pool replicates across dirs (last arg is the output dir):
./1_combine_L001_L002_multiple_replicates.sh <dir1> <dir2> "$BASE_DIR/raw/all_merged_fastqs"

# Stage 2 — QC + multiqc
./2_qc.sh "$BASE_DIR/raw/all_merged_fastqs"

# Stage 3 — host removal
./3_remove_host_reads.sh

# Stage 4 — Kraken2 + Bracken
./4_kraken.sh

# Stage 5 — MetaPhlAn
./5_metaphlan.sh
```

SLURM equivalents:

```bash
# Stage 0
./0_submit_extract_jobs.sh

# Stage 2 — generate the sample list, set the array size it prints, then submit
./make_qc_list.sh "$BASE_DIR/raw/all_merged_fastqs"
#   edit --array=0-(N-1) in 2_qc_array.slurm to match, then:
sbatch 2_qc_array.slurm
#   after the array finishes, build the summary:
multiqc "$BASE_DIR/raw/all_merged_fastqs/fastqc_reports" -o "$BASE_DIR/raw/all_merged_fastqs/multiqc_summary"

# Stage 3 — one job per sample under raw/<subfold>
./3_remove_host_reads_job_submit.sh <subfold>

# Stage 5 — one job per sample
./5_submit_metaphlan_jobs.sh
```

## Notes

- **Idempotent stages.** Stages 2–5 skip a sample when its expected output
  already exists, so re-running after a partial failure only processes what's
  missing.
- **`SUBFOLD` consistency.** Stages 3 and 5 read from and write to
  `<dir>/$SUBFOLD`. Keep `SUBFOLD` set consistently across a run, or set it to
  `""` everywhere if you are not organizing data into batch subdirectories.
- **`#SBATCH` directives are static.** SLURM can't read shell variables in
  `#SBATCH` lines, so the resource values in `2_qc_array.slurm` and
  `3_remove_host_reads_job.sh` are hardcoded there. If you change the
  `SLURM_*` defaults in `config.sh`, update those two files to match.
- **Disk cleanup.** Stages 3, 4, and 5 remove their large intermediates
  (SAM/BAM, decompressed fastq, combined fastq) after each sample.
- **Sample base names** are derived by stripping known suffixes
  (`_R1_combined.fastq.gz`, `_R1_combined_fastp.fastq.gz`,
  `_host_removed_R1.fastq.gz`) — no hand-maintained sample lists are needed
  except for the SLURM array in stage 2, which `make_qc_list.sh` generates.
```

## Tool requirements

`tar`, `fastp`, `fastqc`, `multiqc`, `bowtie2`, `samtools`, `kraken2`,
`bracken`, `metaphlan`, plus a Kraken2/Bracken database, a bowtie2 host index,
and a MetaPhlAn bowtie2 database.
