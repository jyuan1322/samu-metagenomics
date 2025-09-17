#!/bin/bash

# input_dir="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics/raw"
input_dir="$1"
# for filename in $input_dir/SaMu*R1_combined.fastq.gz;do
for filename in $input_dir/*R1_combined.fastq.gz;do
    echo $filename
    base=${filename%_R1_combined.fastq.gz}
    file1=${base}_R1_combined.fastq.gz
    file2=${base}_R2_combined.fastq.gz
    if [ ! -e ${base}"_R1_combined_fastp.fastq.gz" ]; then
        echo $file1
        echo $file2
        fastp \
            --in1 $file1 \
            --out1 "${base}"_R1_combined_fastp.fastq.gz \
            --in2 $file2 \
            --out2 "${base}"_R2_combined_fastp.fastq.gz \
            --detect_adapter_for_pe \
            --thread 8 \
            --html "$base"_fastp.html \
            --json "$base"_fastp.json
    fi
    if [ -e ${base}"_R1_combined_fastp.fastq.gz" ]; then
        if [ ! -e $file1"_fastp_fastqc.html" ]; then
            echo "Checking quality using fastqc for: $filename'_fastp'"
            mkdir -p $input_dir/fastqc_reports
            fastqc \
                --threads 8 \
                --outdir $input_dir/fastqc_reports \
                "${base}"_R*_combined_fastp.fastq.gz
        fi
    fi
done
multiqc $input_dir/fastqc_reports -o $input_dir/multiqc_summary
