#!/bin/bash
# =============================================================================
# 4_kraken.sh — Kraken2 classification + Bracken abundance re-estimation for
# every sample's host-removed reads, with the Bracken output sorted by
# abundance. Plain-loop version (no SLURM submitter in the original set).
#
# Config: BASE_DIR, KRAKEN_DB, CPUS, KRAKEN_*, BRACKEN_*  (see config.sh)
# Usage:  ./4_kraken.sh
#         ./4_kraken.sh 2>&1 | tee logs/4_kraken.log   # to save logs
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

mkdir -p "$KRAKEN_DIR"

for filename in "$HOST_REMOVED_DIR"/*_host_removed_R1.fastq.gz; do
    echo "$filename"
    file=$(basename "$filename")
    base=${file%_host_removed_R1.fastq.gz}

    if [ -f "$KRAKEN_DIR/sorted_${base}_profile.bracken" ]; then
        echo "Skipping because $KRAKEN_DIR/sorted_${base}_profile.bracken exists"
        continue
    fi

    # kraken2 needs decompressed fastq; gunzip -c keeps the .gz originals.
    gunzip -c "$HOST_REMOVED_DIR/${base}_host_removed_R1.fastq.gz" > "$HOST_REMOVED_DIR/${base}_host_removed_R1.fastq"
    gunzip -c "$HOST_REMOVED_DIR/${base}_host_removed_R2.fastq.gz" > "$HOST_REMOVED_DIR/${base}_host_removed_R2.fastq"

    kraken2 \
        --db "$KRAKEN_DB" \
        --paired \
        --threads "$CPUS" \
        --confidence "$KRAKEN_CONFIDENCE" \
        --minimum-base-quality "$KRAKEN_MIN_BASEQ" \
        --report-minimizer-data \
        --report "$KRAKEN_DIR/${base}_profile.report" \
        --output "$KRAKEN_DIR/${base}_kraken_output.txt" \
        "$HOST_REMOVED_DIR/${base}_host_removed_R1.fastq" \
        "$HOST_REMOVED_DIR/${base}_host_removed_R2.fastq"

    bracken -d "$KRAKEN_DB" \
        -i "$KRAKEN_DIR/${base}_profile.report" \
        -o "$KRAKEN_DIR/${base}_profile.bracken" \
        -r "$BRACKEN_READ_LEN" \
        -l "$BRACKEN_LEVEL" \
        -t "$BRACKEN_THRESHOLD"

    # Sort bracken rows by fraction of reads (column 7) descending, keep header.
    (
        head -n 1 "$KRAKEN_DIR/${base}_profile.bracken"
        tail -n +2 "$KRAKEN_DIR/${base}_profile.bracken" | sort -t$'\t' -k7,7nr
    ) > "$KRAKEN_DIR/sorted_${base}_profile.bracken"

    # Remove the decompressed fastq to save space.
    rm -f "$HOST_REMOVED_DIR/${base}_host_removed_R1.fastq" \
          "$HOST_REMOVED_DIR/${base}_host_removed_R2.fastq"
done
