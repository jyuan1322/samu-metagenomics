#!/bin/bash

input_dir="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics"
# subfold="block1" # subfolder under raw, if data is organized in blocks
subfold="$1"

mkdir -p "$input_dir/mapping/${subfold}" "$input_dir/host_removed/${subfold}" "$input_dir/logs"

for r1_file in "$input_dir"/raw/"$subfold"/*_R1_combined_fastp.fastq.gz; do
    base=$(basename "$r1_file" _R1_combined_fastp.fastq.gz)
    echo "Submitting SLURM job for $base"
    sbatch 3_remove_host_reads_job.sh "$base" "$subfold"
done

