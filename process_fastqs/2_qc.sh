#!/bin/bash
# =============================================================================
# 2_qc.sh — adapter/quality trimming (fastp) + quality reports (fastqc) for
# every sample in a directory, then a multiqc summary. Plain-loop version;
# use 2_qc_array.slurm for the SLURM array version.
#
# Config: CPUS  (see config.sh)
# Usage:  ./2_qc.sh <input_dir>
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

input_dir="$1"

for filename in "$input_dir"/*R1_combined.fastq.gz; do
    echo "$filename"
    base=${filename%_R1_combined.fastq.gz}
    file1=${base}_R1_combined.fastq.gz
    file2=${base}_R2_combined.fastq.gz

    if [ ! -e "${base}_R1_combined_fastp.fastq.gz" ]; then
        fastp \
            --in1 "$file1" \
            --out1 "${base}_R1_combined_fastp.fastq.gz" \
            --in2 "$file2" \
            --out2 "${base}_R2_combined_fastp.fastq.gz" \
            --detect_adapter_for_pe \
            --thread "$CPUS" \
            --html "${base}_fastp.html" \
            --json "${base}_fastp.json"
    fi

    if [ -e "${base}_R1_combined_fastp.fastq.gz" ]; then
        if [ ! -e "${base}_R1_combined_fastp_fastqc.html" ]; then
            echo "Running fastqc for: ${base}"
            mkdir -p "$input_dir/fastqc_reports"
            fastqc \
                --threads "$CPUS" \
                --outdir "$input_dir/fastqc_reports" \
                "${base}"_R*_combined_fastp.fastq.gz
        fi
    fi
done

multiqc "$input_dir/fastqc_reports" -o "$input_dir/multiqc_summary"
