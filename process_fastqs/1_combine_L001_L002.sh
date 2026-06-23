#!/bin/bash

# Usage: ./combine_lanes.sh /path/to/input_dir /path/to/output_dir

input_dir="$1"
output_dir="$2"

mkdir -p "$output_dir"

# Get all unique sample prefixes
samples=$(find "$input_dir" -name "*L001_R1_001.fastq.gz" | sed -E 's/(_S[0-9]+)?_L001_R1_001.fastq.gz//' | xargs -n1 basename | sort | uniq)

for sample in $samples; do
    echo "Processing sample: $sample"

    # Construct paths for each lane and read
    r1_l001="${input_dir}/${sample}_S*_L001_R1_001.fastq.gz"
    r1_l002="${input_dir}/${sample}_S*_L002_R1_001.fastq.gz"
    r2_l001="${input_dir}/${sample}_S*_L001_R2_001.fastq.gz"
    r2_l002="${input_dir}/${sample}_S*_L002_R2_001.fastq.gz"

    # Output files
    r1_out="${output_dir}/${sample}_R1_combined.fastq.gz"
    r2_out="${output_dir}/${sample}_R2_combined.fastq.gz"

    # Concatenate matching files for each read direction
    cat $(ls $r1_l001 $r1_l002 2>/dev/null | sort) > "$r1_out"
    cat $(ls $r2_l001 $r2_l002 2>/dev/null | sort) > "$r2_out"

    echo "→ Created: $r1_out and $r2_out"
done

