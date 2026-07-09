#!/bin/bash

# remove the counts column from metaphlan files for compatibility with humann

input_dir="/data/bwh-comppath-full/gerberlab/jy1008/SaMu/metagenomics/metaphlan_v31_out/all_merged_fastqs"
output_dir="/data/bwh-comppath-full/gerberlab/jy1008/SaMu/metagenomics/metaphlan_v31_out/all_merged_fastqs_nocounts"
mkdir -p "$output_dir"

for f in "$input_dir"/*_profile_v31.txt; do
    base=$(basename "$f" _profile_v31.txt)
    out_file="$output_dir/${base}_profile_v31_nocounts.txt"

    awk -F'\t' 'BEGIN{OFS="\t"}
      /^#/ && $0 !~ /^#clade_name/ {print; next}
      /^#clade_name/ {print $1,$2,$3,$4; next}
      {print $1,$2,$3,$4}' "$f" > "$out_file"

    echo "Processed: $base"
done
