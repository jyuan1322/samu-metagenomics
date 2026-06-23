#!/bin/bash
# =============================================================================
# config.sh — central configuration for the metagenomics preprocessing
# pipeline (tar extraction -> lane merging -> QC -> host removal ->
# kraken/bracken -> metaphlan).
#
# Source this from every script:
#   source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
#
# Every value below can also be overridden by exporting an env var of the
# same name before running a script, e.g.:
#   SUBFOLD=block1 bash run/run_2_qc.sh local
# =============================================================================

# ---------------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------------
# BASE_DIR is the project root for this study's data. Everything else is
# derived from it. Point this at your study before running anything.
BASE_DIR="${BASE_DIR:-/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics}"

RAW_DIR="$BASE_DIR/raw"                # combined + QC'd reads
MAPPING_DIR="$BASE_DIR/mapping"        # bowtie2 host-mapping intermediates
HOST_REMOVED_DIR="$BASE_DIR/host_removed"
KRAKEN_DIR="$BASE_DIR/kraken_out"
METAPHLAN_DIR="$BASE_DIR/metaphlan_out"
LOG_DIR="$BASE_DIR/logs"

# SUBFOLD organizes a sequencing batch/run as a subdirectory under raw/,
# mapping/, host_removed/, kraken_out/, metaphlan_out/. Set to "" to operate
# directly on the top-level directories instead.
SUBFOLD="${SUBFOLD:-all_merged_fastqs}"

# ---------------------------------------------------------------------------
# Stage 0: raw tar archives to extract fastq.gz files from
# ---------------------------------------------------------------------------
# Edit this list per study. Each entry is extracted into its own
# RAW_DIR/<tar basename>/ directory by run_0_extract_fastqs.sh.
TAR_FILES=(
  # "/path/to/run1.tar"
  # "/path/to/run2.tar"
)

# ---------------------------------------------------------------------------
# Reference databases
# ---------------------------------------------------------------------------
HOST_BOWTIE2_INDEX="${HOST_BOWTIE2_INDEX:-/data/bwh-comppath-seq/databases/GRCh38_noalt_as/GRCh38_noalt_as}"
KRAKEN_DB="${KRAKEN_DB:-/data/bwh-comppath-seq/databases/uhgg_kraken}"
METAPHLAN_BOWTIE2_DB="${METAPHLAN_BOWTIE2_DB:-/data/bwh-comppath-seq/databases/bowtie2}"

# ---------------------------------------------------------------------------
# Tool parameters
# ---------------------------------------------------------------------------
BRACKEN_READ_LEN="${BRACKEN_READ_LEN:-150}"
BRACKEN_LEVEL="${BRACKEN_LEVEL:-S}"
BRACKEN_THRESHOLD="${BRACKEN_THRESHOLD:-10}"
KRAKEN_CONFIDENCE="${KRAKEN_CONFIDENCE:-0.1}"
KRAKEN_MIN_BASEQ="${KRAKEN_MIN_BASEQ:-20}"

# ---------------------------------------------------------------------------
# Compute resources
# ---------------------------------------------------------------------------
CPUS="${CPUS:-8}"

# SLURM defaults, used by run/*.sh wrappers when invoked in "slurm" mode.
SLURM_PARTITION="${SLURM_PARTITION:-bwh_comppath}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"            # leave empty if not required

SLURM_EXTRACT_TIME="${SLURM_EXTRACT_TIME:-3:00:00}"
SLURM_EXTRACT_MEM="${SLURM_EXTRACT_MEM:-20G}"

SLURM_QC_TIME="${SLURM_QC_TIME:-3:00:00}"
SLURM_QC_MEM="${SLURM_QC_MEM:-16G}"

SLURM_HOST_REMOVAL_TIME="${SLURM_HOST_REMOVAL_TIME:-4:00:00}"
SLURM_HOST_REMOVAL_MEM="${SLURM_HOST_REMOVAL_MEM:-32G}"

SLURM_KRAKEN_TIME="${SLURM_KRAKEN_TIME:-4:00:00}"
SLURM_KRAKEN_MEM="${SLURM_KRAKEN_MEM:-32G}"

SLURM_METAPHLAN_TIME="${SLURM_METAPHLAN_TIME:-8:00:00}"
SLURM_METAPHLAN_MEM="${SLURM_METAPHLAN_MEM:-32G}"

# ---------------------------------------------------------------------------
# Conda / micromamba environments
# ---------------------------------------------------------------------------
MAMBA_EXE="${MAMBA_EXE:-/PHShome/jy1008/bin/micromamba}"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/PHShome/jy1008/.local/share/mamba}"
METAGEN_ENV="${METAGEN_ENV:-metagen-env}"        # fastp, fastqc, bowtie2, samtools, kraken2, bracken
METAPHLAN_ENV="${METAPHLAN_ENV:-metaphlan_env}"  # metaphlan
