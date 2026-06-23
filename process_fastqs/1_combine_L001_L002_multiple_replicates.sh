#!/bin/bash
# =============================================================================
# 1_combine_L001_L002_multiple_replicates.sh — concatenate lanes AND
# replicates into one combined R1/R2 pair per sample, pooling files found
# across any number of input directories. The last argument is always the
# output directory.
#
# Usage: ./1_combine_L001_L002_multiple_replicates.sh <in_dir1> [<in_dir2> ...] <output_dir>
# =============================================================================
set -euo pipefail

output_dir="${@: -1}"
input_dirs=("${@:1:$#-1}")
mkdir -p "$output_dir"

# Sample prefixes across all input dirs (strip optional _S<N> + lane/read suffix).
samples=$(find "${input_dirs[@]}" -name "*L001_R1_001.fastq.gz" \
    | sed -E 's/(_S[0-9]+)?_L00[0-9]_R1_001\.fastq\.gz$//' \
    | xargs -n1 basename | sort -u)

for sample in $samples; do
    echo "Processing sample: $sample"

    r1_files=$(find "${input_dirs[@]}" -name "${sample}_S*_L*_R1_001.fastq.gz" | sort)
    r2_files=$(find "${input_dirs[@]}" -name "${sample}_S*_L*_R2_001.fastq.gz" | sort)

    r1_out="${output_dir}/${sample}_R1_combined.fastq.gz"
    r2_out="${output_dir}/${sample}_R2_combined.fastq.gz"

    if [[ -n "$r1_files" ]]; then
        cat $r1_files > "$r1_out"
    fi
    if [[ -n "$r2_files" ]]; then
        cat $r2_files > "$r2_out"
    fi

    echo "-> Created: $r1_out and $r2_out"
done
