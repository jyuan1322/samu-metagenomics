#!/bin/bash
# =============================================================================
# make_qc_list.sh — regenerate the sample list consumed by 2_qc_array.slurm
# by scanning a directory for *_R1_combined.fastq.gz files. Replaces the
# old hand-maintained 2_qc_files_list.txt.
#
# Usage: ./make_qc_list.sh <input_dir> [output_list]
# Prints the line count so you can set --array=0-(N-1) in 2_qc_array.slurm.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

input_dir="$1"
out_list="${2:-$SCRIPT_DIR/2_qc_files_list.txt}"

find "$input_dir" -maxdepth 1 -type f -name "*_R1_combined.fastq.gz" \
    | sort > "$out_list"

n=$(wc -l < "$out_list")
echo "Wrote $n samples to $out_list"
echo "Set --array=0-$((n-1)) in 2_qc_array.slurm"
