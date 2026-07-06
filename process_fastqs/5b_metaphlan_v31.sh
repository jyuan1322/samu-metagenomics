#!/bin/bash
# =============================================================================
# 5b_metaphlan_v31.sh — MetaPhlAn profiling against the pre-SGB v31 database
# (mpa_v31_CHOCOPhlAn_201901, MetaPhlAn 3.1.0), for every sample's
# host-removed reads. This profile is what HUMAnN's --taxonomic-profile
# consumes in 6_humann.sh — it must NOT be confused with the vJun23 profiles
# from 5_metaphlan.sh, which use SGB-based naming incompatible with HUMAnN's
# ChocoPhlAn.
# Plain-loop version; use 5b_submit_metaphlan_v31_jobs.sh for the SLURM version.
#
# Config: BASE_DIR, SUBFOLD, METAPHLAN_V31_*, CPUS  (see config.sh)
# Usage:  ./5b_metaphlan_v31.sh
#         ./5b_metaphlan_v31.sh 2>&1 | tee logs/5b_metaphlan_v31.log
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

host_removed_sub="$HOST_REMOVED_DIR/$SUBFOLD"
output_dir="$METAPHLAN_V31_DIR/$SUBFOLD"
mkdir -p "$output_dir"

for r1_file in "$host_removed_sub"/*_host_removed_R1.fastq.gz; do
    base=$(basename "$r1_file" _host_removed_R1.fastq.gz)
    r2_file="$host_removed_sub/${base}_host_removed_R2.fastq.gz"
    comb_file="$host_removed_sub/${base}_host_removed_R1R2_combined_v31.fastq.gz"

    out_file="$output_dir/${base}_profile_v31.txt"
    if [[ -f "$out_file" ]]; then
        echo "Skipping $base (already completed)"
        continue
    fi

    echo "Combining reads for $base"
    gunzip -c "$r1_file" "$r2_file" | gzip > "$comb_file"

    echo "Running MetaPhlAn (v31, pre-SGB) on $base"
    metaphlan \
        "$comb_file" \
        --input_type fastq \
        --bowtie2db "$METAPHLAN_V31_BOWTIE2_DB" \
        --index "$METAPHLAN_V31_INDEX" \
        --nproc "$CPUS" \
        -t rel_ab_w_read_stats \
        -o "$out_file"

    echo "Cleaning up combined file"
    rm -f "$comb_file"
done