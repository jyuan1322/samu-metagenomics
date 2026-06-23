#!/bin/bash
# =============================================================================
# 0_collect_tar_fastqs.sh — extract every tar archive in TAR_FILES and collect
# the fastq.gz files into RAW_DIR. Plain-loop version (runs serially in this
# shell); use 0_submit_extract_jobs.sh for the SLURM version.
#
# Config: TAR_FILES, RAW_DIR  (see config.sh)
# Usage:  ./0_collect_tar_fastqs.sh
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
    log "Processing $TAR_FILE..."

    if [[ ! -f "$TAR_FILE" ]]; then
        echo "  [Warning] File not found: $TAR_FILE" >&2
        continue
    fi

    tar -xvf "$TAR_FILE"

    # Move any extracted fastq.gz (up to 2 levels deep) into RAW_DIR.
    find . -maxdepth 2 -type f -name "*.fastq.gz" -exec mv -v {} "$RAW_DIR/" \;
done

log "Done. All .fastq.gz files have been moved to $RAW_DIR/"
