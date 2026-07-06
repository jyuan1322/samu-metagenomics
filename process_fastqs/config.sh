#!/bin/bash
# =============================================================================
# config_test.sh — test configuration for a single-sample dry run of the
# MetaPhlAn v31 -> HUMAnN pipeline. Mirrors config.sh, with BASE_DIR and
# SUBFOLD pointed at a small manually-created test directory instead of the
# real study data.
#
# Usage: source this INSTEAD of config.sh, e.g.:
#   CONFIG=config_test.sh bash 5b_metaphlan_v31.sh
# (only works if your scripts source $CONFIG rather than a hardcoded
#  config.sh — see note at bottom if they don't yet)
# =============================================================================

# ---------------------------------------------------------------------------
# Directory layout — test-specific
# ---------------------------------------------------------------------------
BASE_DIR="/data/bwh-comppath-full/databases/test_run"

RAW_DIR="$BASE_DIR/raw"
MAPPING_DIR="$BASE_DIR/mapping"
HOST_REMOVED_DIR="/data/bwh-comppath-full/databases/test_SaMu9_host_removed"
KRAKEN_DIR="$BASE_DIR/kraken_out"
METAPHLAN_DIR="$BASE_DIR/metaphlan_out"
LOG_DIR="$BASE_DIR/logs"

# No subfolder structure for this test — files sit directly in HOST_REMOVED_DIR
SUBFOLD=""

TAR_FILES=()

# ---------------------------------------------------------------------------
# Reference databases — same as production, no need to duplicate downloads
# ---------------------------------------------------------------------------
HOST_BOWTIE2_INDEX="/data/bwh-comppath-seq/databases/GRCh38_noalt_as/GRCh38_noalt_as"
KRAKEN_DB="/data/bwh-comppath-seq/databases/uhgg_kraken"
METAPHLAN_BOWTIE2_DB="/data/bwh-comppath-seq/databases/bowtie2"

# ---------------------------------------------------------------------------
# Tool parameters — unchanged from production
# ---------------------------------------------------------------------------
BRACKEN_READ_LEN="150"
BRACKEN_LEVEL="S"
BRACKEN_THRESHOLD="10"
KRAKEN_CONFIDENCE="0.1"
KRAKEN_MIN_BASEQ="20"

# ---------------------------------------------------------------------------
# Compute resources — lighter, since this is a tiny subsampled test
# ---------------------------------------------------------------------------
CPUS="2"

SLURM_PARTITION="bwh_comppath"
SLURM_ACCOUNT=""

SLURM_METAPHLAN_TIME="0:30:00"
SLURM_METAPHLAN_MEM="8G"

SLURM_METAPHLAN_V31_TIME="0:30:00"
SLURM_METAPHLAN_V31_MEM="8G"

SLURM_HUMANN_TIME="1:00:00"
SLURM_HUMANN_MEM="16G"

# ---------------------------------------------------------------------------
# Conda / micromamba environments — same as production
# ---------------------------------------------------------------------------
MAMBA_EXE="/PHShome/jy1008/bin/micromamba"
MAMBA_ROOT_PREFIX="/PHShome/jy1008/.local/share/mamba"
METAGEN_ENV="metagen-env"
METAPHLAN_ENV="metaphlan_env"

# ---------------------------------------------------------------------------
# MetaPhlAn v31 (pre-SGB) + HUMAnN — same databases/envs as production
# ---------------------------------------------------------------------------
METAPHLAN_V31_ENV="metaphlan31_env"
METAPHLAN_V31_BOWTIE2_DB="/data/bwh-comppath-full/databases/humann/metaphlan_v31"
METAPHLAN_V31_INDEX="mpa_v31_CHOCOPhlAn_201901"
METAPHLAN_V31_DIR="$BASE_DIR/metaphlan_v31_out"

HUMANN_DIR="$BASE_DIR/humann_out"
HUMANN_CHOCOPHLAN_DB="/data/bwh-comppath-full/databases/humann/chocophlan/chocophlan"
HUMANN_UNIREF_DB="/data/bwh-comppath-full/databases/humann/uniref/uniref"
HUMANN_ENV="humann39_env"