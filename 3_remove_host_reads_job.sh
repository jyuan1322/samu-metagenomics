#!/bin/bash
#SBATCH --job-name=host_removal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/data/bwh-comppath-seq/jy1008/SaMu/logs/%x_%j.out
#SBATCH --error=/data/bwh-comppath-seq/jy1008/SaMu/logs/%x_%j.err

# Initialize micromamba for this shell
eval "$(micromamba shell hook --shell bash)"
micromamba activate metagen-env

set -euo pipefail

input_dir="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics"
cpus=8

base="$1"
subfold="$2"

r1_file="$input_dir/raw/${subfold}/${base}_R1_combined_fastp.fastq.gz"
r2_file="$input_dir/raw/${subfold}/${base}_R2_combined_fastp.fastq.gz"

out_bam="$input_dir/mapping/${subfold}/${base}_bothReadsUnmapped_sorted.bam"

if [[ -e "$r1_file" && -e "$r2_file" && ! -e "$out_bam" ]]; then
    echo "Running host removal for sample: $base"

    bowtie2 -p $cpus -x /data/bwh-comppath-seq/databases/GRCh38_noalt_as/GRCh38_noalt_as \
        -1 "$r1_file" -2 "$r2_file" \
        -S "$input_dir/mapping/${subfold}/${base}_mapped_and_unmapped.sam"

    samtools view -bS "$input_dir/mapping/${subfold}/${base}_mapped_and_unmapped.sam" \
        > "$input_dir/mapping/${subfold}/${base}_mapped_and_unmapped.bam"

    samtools view -b -f 12 -F 256 "$input_dir/mapping/${subfold}/${base}_mapped_and_unmapped.bam" \
        > "$input_dir/mapping/${subfold}/${base}_bothReadsUnmapped.bam"

    samtools sort -n -m 5G -@ $cpus "$input_dir/mapping/${subfold}/${base}_bothReadsUnmapped.bam" \
        -o "$out_bam"

    samtools fastq -@ $cpus "$out_bam" \
        -1 "$input_dir/host_removed/${subfold}/${base}_host_removed_R1.fastq.gz" \
        -2 "$input_dir/host_removed/${subfold}/${base}_host_removed_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n

    # Optional cleanup
    rm "$input_dir/mapping/${subfold}/${base}_mapped_and_unmapped.sam" \
       "$input_dir/mapping/${subfold}/${base}_mapped_and_unmapped.bam" \
       "$input_dir/mapping/${subfold}/${base}_bothReadsUnmapped.bam"

else
    echo "Input files missing or output already exists for $base, skipping."
fi

