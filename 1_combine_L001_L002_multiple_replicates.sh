#!/bin/bash

# Usage: ./1_combine_L001_L002_multiple_replicates.sh /path/to/input_dir1 /path/to/input_dir2 ... /path/to/output_dir
# The last argument is always the output directory.

# Extract output dir (last argument)
output_dir="${@: -1}"
mkdir -p "$output_dir"

# All input dirs are everything except last argument
input_dirs=("${@:1:$#-1}")

# Find all sample prefixes across all input dirs
samples=$(find "${input_dirs[@]}" -name "*L001_R1_001.fastq.gz" \
    | sed -E 's/(_S[0-9]+)?_L00[0-9]_R1_001.fastq.gz//' \
    | xargs -n1 basename | sort | uniq)

for sample in $samples; do
    echo "Processing sample: $sample"

    # Collect all files for this sample across all dirs
    r1_files=$(find "${input_dirs[@]}" -name "${sample}_S*_L*_R1_001.fastq.gz" | sort)
    r2_files=$(find "${input_dirs[@]}" -name "${sample}_S*_L*_R2_001.fastq.gz" | sort)

    # Output files
    r1_out="${output_dir}/${sample}_R1_combined.fastq.gz"
    r2_out="${output_dir}/${sample}_R2_combined.fastq.gz"

    # Concatenate across replicates and lanes
    if [[ -n "$r1_files" ]]; then
        cat $r1_files > "$r1_out"
    fi
    if [[ -n "$r2_files" ]]; then
        cat $r2_files > "$r2_out"
    fi

    echo "→ Created: $r1_out and $r2_out"
done

