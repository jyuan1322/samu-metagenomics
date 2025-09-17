#!/bin/bash
#SBATCH --job-name=fastp_qc
#SBATCH --output=logs/fastp_qc_%A_%a.out
#SBATCH --error=logs/fastp_qc_%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-201   # adjust based on number of samples in samples.txt

# Initialize micromamba for this shell
eval "$(micromamba shell hook --shell bash)"
micromamba activate metagen-env

set -euo pipefail

# Input setup
input_dir="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics/raw/all_merged_fastqs"
sample_list="/data/bwh-comppath-seq/jy1008/SaMu/2_qc_files_list.txt"

# Pick the sample for this array task
filename=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $sample_list)
base=${filename%_R1_combined.fastq.gz}
file1=${base}_R1_combined.fastq.gz
file2=${base}_R2_combined.fastq.gz

# Run fastp
if [ ! -e ${base}"_R1_combined_fastp.fastq.gz" ]; then
    echo "Running fastp for $base"
    fastp \
        --in1 $file1 \
        --out1 "${base}"_R1_combined_fastp.fastq.gz \
        --in2 $file2 \
        --out2 "${base}"_R2_combined_fastp.fastq.gz \
        --detect_adapter_for_pe \
        --thread $SLURM_CPUS_PER_TASK \
        --html "$base"_fastp.html \
        --json "$base"_fastp.json
fi

# Run fastqc
if [ -e ${base}"_R1_combined_fastp.fastq.gz" ]; then
    mkdir -p $input_dir/fastqc_reports
    fastqc \
        --threads $SLURM_CPUS_PER_TASK \
        --outdir $input_dir/fastqc_reports \
        "${base}"_R*_combined_fastp.fastq.gz
fi

