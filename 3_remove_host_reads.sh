#!/bin/bash

input_dir="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics"

cpus=8
mkdir -p "$input_dir/mapping" "$input_dir/host_removed"

for filename in "$input_dir"/raw/SaMu*_R1_combined_fastp.fastq.gz; do
    echo "$filename"
    folder=$(dirname "$filename")
    file=$(basename "$filename")
    base=${file%_R1_combined_fastp.fastq.gz}
    file2=${base}_R2_combined_fastp.fastq.gz

    if [ ! -e "$input_dir/mapping/${base}_bothReadsUnmapped_sorted.bam" ]; then
        echo "Aligning ${base} with bowtie2"
        bowtie2 -p $cpus -x /data/bwh-comppath-seq/databases/GRCh38_noalt_as/GRCh38_noalt_as \
            -1 "$folder/$file" -2 "$folder/$file2" \
            -S "$input_dir/mapping/${base}_mapped_and_unmapped.sam"

        samtools view -bS "$input_dir/mapping/${base}_mapped_and_unmapped.sam" > "$input_dir/mapping/${base}_mapped_and_unmapped.bam"
        samtools view -b -f 12 -F 256 "$input_dir/mapping/${base}_mapped_and_unmapped.bam" > "$input_dir/mapping/${base}_bothReadsUnmapped.bam"
        samtools sort -n -m 5G -@ $cpus "$input_dir/mapping/${base}_bothReadsUnmapped.bam" -o "$input_dir/mapping/${base}_bothReadsUnmapped_sorted.bam"

        samtools fastq -@ $cpus "$input_dir/mapping/${base}_bothReadsUnmapped_sorted.bam" \
            -1 "$input_dir/host_removed/${base}_host_removed_R1.fastq.gz" \
            -2 "$input_dir/host_removed/${base}_host_removed_R2.fastq.gz" \
            -0 /dev/null -s /dev/null -n
    fi
done
