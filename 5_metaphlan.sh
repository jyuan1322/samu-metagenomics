#!/bin/bash
set -euo pipefail

# run like: ./4_metaphlan.sh 2>&1 | tee logs/4_metaphlan.log

cpus=8
input_dir="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics"
subfold="all_merged_fastqs"
output_dir="$input_dir/metaphlan_out/$subfold"
mkdir -p "$output_dir"

for r1_file in "$input_dir"/host_removed/"$subfold"/*_host_removed_R1.fastq.gz; do
    base=$(basename "$r1_file" _host_removed_R1.fastq.gz)
    r2_file="$input_dir/host_removed/$subfold/${base}_host_removed_R2.fastq.gz"
    comb_file="$input_dir/host_removed/$subfold/${base}_host_removed_R1R2_combined.fastq.gz"

    echo "Combining reads for $base"
    gunzip -c "$r1_file" "$r2_file" | gzip > "$comb_file"

    echo "Running MetaPhlAn on $base"
    metaphlan \
        "$comb_file" \
        --input_type fastq \
        --bowtie2db /data/bwh-comppath-seq/databases/bowtie2/ \
        --nproc "$cpus" \
        -o "$output_dir/${base}_profile.txt"

    echo "Cleaning up combined file"
    rm "$comb_file"
done
