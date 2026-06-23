#!/bin/bash
#SBATCH --job-name=host_removal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# =============================================================================
# 3_remove_host_reads_job.sh — host removal for ONE sample. Submitted once per
# sample by 3_remove_host_reads_job_submit.sh. Data is organized under an
# optional SUBFOLD layer (e.g. a sequencing batch) within raw/mapping/
# host_removed.
#
# NOTE: #SBATCH directives can't read shell variables; keep them in sync with
# config.sh if you change defaults. The --output/--error paths are relative,
# so submit from the pipeline directory (which should contain a logs/ dir).
#
# Usage (via sbatch): 3_remove_host_reads_job.sh <sample_base> <subfold>
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

activate_env "$METAGEN_ENV"

base="$1"
subfold="$2"

raw_sub="$RAW_DIR/${subfold}"
mapping_sub="$MAPPING_DIR/${subfold}"
host_removed_sub="$HOST_REMOVED_DIR/${subfold}"
mkdir -p "$mapping_sub" "$host_removed_sub"

r1_file="$raw_sub/${base}_R1_combined_fastp.fastq.gz"
r2_file="$raw_sub/${base}_R2_combined_fastp.fastq.gz"
out_bam="$mapping_sub/${base}_bothReadsUnmapped_sorted.bam"

if [[ -e "$r1_file" && -e "$r2_file" && ! -e "$out_bam" ]]; then
    echo "Running host removal for sample: $base"

    bowtie2 -p "$CPUS" -x "$HOST_BOWTIE2_INDEX" \
        -1 "$r1_file" -2 "$r2_file" \
        -S "$mapping_sub/${base}_mapped_and_unmapped.sam"

    samtools view -bS "$mapping_sub/${base}_mapped_and_unmapped.sam" \
        > "$mapping_sub/${base}_mapped_and_unmapped.bam"
    samtools view -b -f 12 -F 256 "$mapping_sub/${base}_mapped_and_unmapped.bam" \
        > "$mapping_sub/${base}_bothReadsUnmapped.bam"
    samtools sort -n -m 5G -@ "$CPUS" "$mapping_sub/${base}_bothReadsUnmapped.bam" \
        -o "$out_bam"

    samtools fastq -@ "$CPUS" "$out_bam" \
        -1 "$host_removed_sub/${base}_host_removed_R1.fastq.gz" \
        -2 "$host_removed_sub/${base}_host_removed_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n

    rm -f "$mapping_sub/${base}_mapped_and_unmapped.sam" \
          "$mapping_sub/${base}_mapped_and_unmapped.bam" \
          "$mapping_sub/${base}_bothReadsUnmapped.bam"
else
    echo "Input files missing or output already exists for $base, skipping."
fi
