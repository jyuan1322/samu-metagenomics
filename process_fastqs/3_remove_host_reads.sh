#!/bin/bash
# =============================================================================
# 3_remove_host_reads.sh — align each sample to the host genome (bowtie2) and
# keep read pairs where both mates are unmapped (microbial reads). Plain-loop
# version; use 3_remove_host_reads_job_submit.sh for the SLURM version.
#
# Config: BASE_DIR, HOST_BOWTIE2_INDEX, CPUS  (see config.sh)
# Usage:  ./3_remove_host_reads.sh
# Operates on $BASE_DIR/raw, writing to $BASE_DIR/mapping and
# $BASE_DIR/host_removed.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

mkdir -p "$MAPPING_DIR" "$HOST_REMOVED_DIR"

for filename in "$RAW_DIR"/*_R1_combined_fastp.fastq.gz; do
    echo "$filename"
    folder=$(dirname "$filename")
    file=$(basename "$filename")
    base=${file%_R1_combined_fastp.fastq.gz}
    file2=${base}_R2_combined_fastp.fastq.gz

    out_bam="$MAPPING_DIR/${base}_bothReadsUnmapped_sorted.bam"
    if [ -e "$out_bam" ]; then
        continue
    fi

    echo "Aligning ${base} with bowtie2"
    bowtie2 -p "$CPUS" -x "$HOST_BOWTIE2_INDEX" \
        -1 "$folder/$file" -2 "$folder/$file2" \
        -S "$MAPPING_DIR/${base}_mapped_and_unmapped.sam"

    samtools view -bS "$MAPPING_DIR/${base}_mapped_and_unmapped.sam" \
        > "$MAPPING_DIR/${base}_mapped_and_unmapped.bam"
    # -f 12: both reads unmapped. -F 256: exclude secondary alignments.
    samtools view -b -f 12 -F 256 "$MAPPING_DIR/${base}_mapped_and_unmapped.bam" \
        > "$MAPPING_DIR/${base}_bothReadsUnmapped.bam"
    samtools sort -n -m 5G -@ "$CPUS" "$MAPPING_DIR/${base}_bothReadsUnmapped.bam" \
        -o "$out_bam"

    samtools fastq -@ "$CPUS" "$out_bam" \
        -1 "$HOST_REMOVED_DIR/${base}_host_removed_R1.fastq.gz" \
        -2 "$HOST_REMOVED_DIR/${base}_host_removed_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n
done
