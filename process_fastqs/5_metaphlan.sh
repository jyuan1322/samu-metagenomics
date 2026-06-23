#!/bin/bash
# =============================================================================
# 5_metaphlan.sh — MetaPhlAn profiling for every sample's host-removed reads.
# Plain-loop version; use 5_submit_metaphlan_jobs.sh for the SLURM version.
#
# Config: BASE_DIR, SUBFOLD, METAPHLAN_BOWTIE2_DB, CPUS  (see config.sh)
# Usage:  ./5_metaphlan.sh
#         ./5_metaphlan.sh 2>&1 | tee logs/5_metaphlan.log
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

host_removed_sub="$HOST_REMOVED_DIR/$SUBFOLD"
output_dir="$METAPHLAN_DIR/$SUBFOLD"
mkdir -p "$output_dir"

for r1_file in "$host_removed_sub"/*_host_removed_R1.fastq.gz; do
    base=$(basename "$r1_file" _host_removed_R1.fastq.gz)
    r2_file="$host_removed_sub/${base}_host_removed_R2.fastq.gz"
    comb_file="$host_removed_sub/${base}_host_removed_R1R2_combined.fastq.gz"

    out_file="$output_dir/${base}_profile.txt"
    if [[ -f "$out_file" ]]; then
        echo "Skipping $base (already completed)"
        continue
    fi

    echo "Combining reads for $base"
    gunzip -c "$r1_file" "$r2_file" | gzip > "$comb_file"

    echo "Running MetaPhlAn on $base"
    metaphlan \
        "$comb_file" \
        --input_type fastq \
        --bowtie2db "$METAPHLAN_BOWTIE2_DB" \
        --nproc "$CPUS" \
        -t rel_ab_w_read_stats \
        -o "$out_file"

    echo "Cleaning up combined file"
    rm -f "$comb_file"
done
