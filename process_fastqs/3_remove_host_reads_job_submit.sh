#!/bin/bash
# =============================================================================
# 3_remove_host_reads_job_submit.sh — submit one host-removal SLURM job per
# sample found under $RAW_DIR/<subfold>. SLURM version of
# 3_remove_host_reads.sh.
#
# Config: BASE_DIR (and derived dirs)  (see config.sh)
# Usage:  ./3_remove_host_reads_job_submit.sh <subfold>
#         <subfold> is the subdirectory under raw/ to process (e.g. a
#         sequencing batch). Pass "" to operate on the top-level raw/ dir.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

subfold="${1:-}"

mkdir -p "$MAPPING_DIR/${subfold}" "$HOST_REMOVED_DIR/${subfold}" "$LOG_DIR"

for r1_file in "$RAW_DIR/$subfold"/*_R1_combined_fastp.fastq.gz; do
    [ -e "$r1_file" ] || { echo "No fastp reads found in $RAW_DIR/$subfold"; exit 1; }
    base=$(basename "$r1_file" _R1_combined_fastp.fastq.gz)
    echo "Submitting SLURM job for $base"
    sbatch "$SCRIPT_DIR/3_remove_host_reads_job.sh" "$base" "$subfold"
done
