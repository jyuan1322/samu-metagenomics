#!/bin/bash
set -euo pipefail

# run the following to view and save logs from stdout
# ./4_kraken.sh 2>&1 | tee logs/4_kraken.log

cpus=8
input_dir="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics"

for filename in $input_dir/host_removed/SaMu*_host_removed_R1.fastq.gz; do
    echo "$filename"

    # filename="$input_dir/host_removed/SaMu11_host_removed_R1.fastq.gz"
    file=$(basename ${filename})
    base=${file%_host_removed_R1.fastq.gz}
    mkdir -p "$input_dir/kraken_out"
    if [ -f "$input_dir/kraken_out/sorted_${base}_profile.bracken" ]; then
        echo "Skipping because $input_dir/kraken_out/sorted_${base}_profile.bracken exists"
        continue
    fi

    gunzip -c "$input_dir/host_removed/${base}_host_removed_R1.fastq.gz" > "$input_dir/host_removed/${base}_host_removed_R1.fastq"
    gunzip -c "$input_dir/host_removed/${base}_host_removed_R2.fastq.gz" > "$input_dir/host_removed/${base}_host_removed_R2.fastq"

    # cat "$input_dir/host_removed/${base}_host_removed_R1.fastq" \
    #     "$input_dir/host_removed/${base}_host_removed_R2.fastq" > \
    #     "$input_dir/host_removed/${base}_host_removed_R1R2_combined.fastq"

    # metaphlan \
    #     --input_type fastq \
    #     "$input_dir/host_removed/${base}_host_removed_R1R2_combined.fastq" \
    #     --bowtie2db /data/bwh-comppath-seq/databases/bowtie2/ > \
    #     metaphlan/"$beginning"_profile.txt

    kraken2 \
      --db /data/bwh-comppath-seq/databases/uhgg_kraken/ \
      --paired \
      --threads $cpus \
      --confidence 0.1 \
      --minimum-base-quality 20 \
      --report-minimizer-data \
      --report "$input_dir/kraken_out/${base}_profile.report" \
      --output "$input_dir/kraken_out/${base}_kraken_output.txt" \
      "$input_dir/host_removed/${base}_host_removed_R1.fastq" \
      "$input_dir/host_removed/${base}_host_removed_R2.fastq"

      # --unclassified-out ${base}_kraken_unclassified.fastq \
      # --classified-out "$base"_kraken_classified.fastq \
      # "$input_dir/host_removed/${base}_host_removed_R1.fastq.gz" \
      # "$input_dir/host_removed/${base}_host_removed_R2.fastq.gz"

    bracken -d /data/bwh-comppath-seq/databases/uhgg_kraken/ \
           -i "$input_dir/kraken_out/${base}_profile.report" \
           -o "$input_dir/kraken_out/${base}_profile.bracken" \
           -r 150 \
           -l S \
           -t 10

    (
        head -n 1 "$input_dir/kraken_out/${base}_profile.bracken" && \
        tail -n +2 "$input_dir/kraken_out/${base}_profile.bracken" | \
        sort -t$'\t' -k7,7nr
    ) > "$input_dir/kraken_out/sorted_${base}_profile.bracken"

done
