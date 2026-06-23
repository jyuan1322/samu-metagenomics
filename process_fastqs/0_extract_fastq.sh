#!/bin/bash
# =============================================================================
# 0_extract_fastq.sh — extract a single tar archive and flatten its fastq.gz
# files into the output directory. Called once per archive by
# 0_submit_extract_jobs.sh (one SLURM job each).
#
# Usage: ./0_extract_fastq.sh <tar_file> <outdir>
# =============================================================================
set -euo pipefail

TAR_FILE="$1"
OUTDIR="$2"

echo "Extracting: $TAR_FILE into $OUTDIR"

mkdir -p "$OUTDIR"

# Extract (preserving the original tar), then flatten fastq.gz into OUTDIR.
tar -xvf "$TAR_FILE" -C "$OUTDIR"
find "$OUTDIR" -mindepth 2 -type f -name "*.fastq.gz" -exec mv -v {} "$OUTDIR" \;

echo "Finished processing $TAR_FILE"
