#!/bin/bash
# =============================================================================
# 6_humann.sh — HUMAnN functional profiling for every sample, using the
# standalone MetaPhlAn v31 (pre-SGB) profile from 5b_metaphlan_v31.sh via
# --taxonomic-profile. No internal MetaPhlAn invocation — HUMAnN's ChocoPhlAn
# (v201901_v31) requires pre-SGB species naming, which only the v31 profile
# provides (NOT the vJun23 profile from 5_metaphlan.sh).
# Plain-loop version; use 6_submit_humann_jobs.sh for the SLURM version.
#
# Config: BASE_DIR, SUBFOLD, HUMANN_*, METAPHLAN_V31_DIR, CPUS  (see config.sh)
# Usage:  ./6_humann.sh
#         ./6_humann.sh 2>&1 | tee logs/6_humann.log
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

host_removed_sub="$HOST_REMOVED_DIR/$SUBFOLD"
profile_sub="$METAPHLAN_V31_DIR/$SUBFOLD"
output_dir="$HUMANN_DIR/$SUBFOLD"
mkdir -p "$output_dir"

for r1_file in "$host_removed_sub"/*_host_removed_R1.fastq.gz; do
    base=$(basename "$r1_file" _host_removed_R1.fastq.gz)
    r2_file="$host_removed_sub/${base}_host_removed_R2.fastq.gz"
    comb_file="$host_removed_sub/${base}_host_removed_R1R2_combined_humann.fastq.gz"
    profile_file="$profile_sub/${base}_profile_v31.txt"

    sample_out_dir="$output_dir/$base"

    if [[ -d "$sample_out_dir" ]] && ls "$sample_out_dir"/*_genefamilies.tsv >/dev/null 2>&1; then
        echo "Skipping $base (already completed)"
        continue
    fi

    if [[ ! -f "$profile_file" ]]; then
        echo "No v31 MetaPhlAn profile found for $base at $profile_file — run 5b_metaphlan_v31.sh first. Skipping."
        continue
    fi

    mkdir -p "$sample_out_dir"

    echo "Combining reads for $base"
    gunzip -c "$r1_file" "$r2_file" | gzip > "$comb_file"

    echo "Running HUMAnN on $base (taxonomic-profile: $profile_file)"
    humann \
        --input "$comb_file" \
        --output "$sample_out_dir" \
        --taxonomic-profile "$profile_file" \
        --nucleotide-database "$HUMANN_CHOCOPHLAN_DB" \
        --protein-database "$HUMANN_UNIREF_DB" \
        --threads "$CPUS"

    echo "Cleaning up combined file"
    rm -f "$comb_file"
done