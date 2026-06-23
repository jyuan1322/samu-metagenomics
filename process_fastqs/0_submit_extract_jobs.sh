#!/bin/bash
# =============================================================================
# 0_submit_extract_jobs.sh — submit one SLURM job per tar archive in
# TAR_FILES; each extracts into its own subdirectory under RAW_DIR. SLURM
# version of 0_collect_tar_fastqs.sh (which it calls via 0_extract_fastq.sh).
#
# Config: TAR_FILES, RAW_DIR, SLURM_PARTITION, SLURM_EXTRACT_TIME/MEM
# Usage:  ./0_submit_extract_jobs.sh
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

mkdir -p "$RAW_DIR"

if [ "${#TAR_FILES[@]}" -eq 0 ]; then
    echo "No TAR_FILES configured in config.sh — add tar archive paths and re-run." >&2
    exit 1
fi

for TAR_FILE in "${TAR_FILES[@]}"; do
    BASENAME=$(basename "$TAR_FILE" .tar)
    OUTDIR="$RAW_DIR/$BASENAME"
    mkdir -p "$OUTDIR"

    log "Submitting extract job for $BASENAME"
    sbatch --job-name="untar_${BASENAME}" \
           --output="$OUTDIR/slurm-%j.out" \
           --error="$OUTDIR/slurm-%j.err" \
           --partition="$SLURM_PARTITION" \
           --time="$SLURM_EXTRACT_TIME" \
           --mem="$SLURM_EXTRACT_MEM" \
           --wrap="bash '$SCRIPT_DIR/0_extract_fastq.sh' '$TAR_FILE' '$OUTDIR'"
done

log "All SLURM jobs submitted."
