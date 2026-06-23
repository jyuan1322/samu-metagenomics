#!/bin/bash
# =============================================================================
# 1_combine_L001_L002.sh — concatenate L001/L002 lanes into one combined
# R1/R2 pair per sample, from a single input directory.
# For data split across multiple replicate directories, use
# 1_combine_L001_L002_multiple_replicates.sh instead.
#
# Usage: ./1_combine_L001_L002.sh <input_dir> <output_dir>
# =============================================================================
set -euo pipefail

input_dir="$1"
output_dir="$2"

mkdir -p "$output_dir"

# Unique sample prefixes (strip optional _S<N> index + lane/read suffix).
samples=$(find "$input_dir" -name "*L001_R1_001.fastq.gz" \
    | sed -E 's/(_S[0-9]+)?_L001_R1_001\.fastq\.gz$//' \
    | xargs -n1 basename | sort -u)

for sample in $samples; do
    echo "Processing sample: $sample"

    r1_l001="${input_dir}/${sample}_S*_L001_R1_001.fastq.gz"
    r1_l002="${input_dir}/${sample}_S*_L002_R1_001.fastq.gz"
    r2_l001="${input_dir}/${sample}_S*_L001_R2_001.fastq.gz"
    r2_l002="${input_dir}/${sample}_S*_L002_R2_001.fastq.gz"

    r1_out="${output_dir}/${sample}_R1_combined.fastq.gz"
    r2_out="${output_dir}/${sample}_R2_combined.fastq.gz"

    cat $(ls $r1_l001 $r1_l002 2>/dev/null | sort) > "$r1_out"
    cat $(ls $r2_l001 $r2_l002 2>/dev/null | sort) > "$r2_out"

    echo "-> Created: $r1_out and $r2_out"
done
