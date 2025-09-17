#!/bin/bash
set -euo pipefail

cpus=8
input_dir="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics"
subfold="all_merged_fastqs"
output_dir="$input_dir/metaphlan_out/$subfold"
log_dir="$input_dir/logs/metaphlan"
mkdir -p "$output_dir" "$log_dir"

for r1_file in "$input_dir"/host_removed/"$subfold"/*_host_removed_R1.fastq.gz; do
    base=$(basename "$r1_file" _host_removed_R1.fastq.gz)
    out_file="$output_dir/${base}_profile.txt"

    # Skip if output already exists
    if [[ -f "$out_file" ]]; then
        echo "Skipping $base (already completed)"
        continue
    fi

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=mp_${base}
#SBATCH --cpus-per-task=$cpus
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=$log_dir/${base}.out
#SBATCH --error=$log_dir/${base}.err

export MAMBA_EXE="/PHShome/jy1008/bin/micromamba"
export MAMBA_ROOT_PREFIX="/PHShome/jy1008/.local/share/mamba"

__mamba_setup="\$("\$MAMBA_EXE" shell hook --shell bash --root-prefix "\$MAMBA_ROOT_PREFIX" 2>/dev/null || true)"
if [ -n "\$__mamba_setup" ]; then
    eval "\$__mamba_setup"
    unset __mamba_setup
else
    echo "Failed to init micromamba"
    exit 1
fi

micromamba activate metaphlan_env || { echo "Failed to activate env"; exit 1; }

echo "Running MetaPhlAn for \$base"

r1_file="$input_dir/host_removed/$subfold/${base}_host_removed_R1.fastq.gz"
r2_file="$input_dir/host_removed/$subfold/${base}_host_removed_R2.fastq.gz"
comb_file="$input_dir/host_removed/$subfold/${base}_host_removed_R1R2_combined.fastq.gz"

gunzip -c "\$r1_file" "\$r2_file" | gzip > "\$comb_file"

metaphlan \\
    "\$comb_file" \\
    --input_type fastq \\
    --bowtie2db /data/bwh-comppath-seq/databases/bowtie2/ \\
    --nproc $cpus \\
    -t rel_ab_w_read_stats \\
    -o "$out_file"

rm "\$comb_file"

echo "Done with \$base"
EOF

done
