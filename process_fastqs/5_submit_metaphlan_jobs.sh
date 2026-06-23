#!/bin/bash
# =============================================================================
# 5_submit_metaphlan_jobs.sh — submit one MetaPhlAn SLURM job per sample under
# $HOST_REMOVED_DIR/$SUBFOLD. SLURM version of 5_metaphlan.sh.
#
# Config: BASE_DIR, SUBFOLD, METAPHLAN_BOWTIE2_DB, CPUS, METAPHLAN_ENV,
#         SLURM_METAPHLAN_*  (see config.sh)
# Usage:  ./5_submit_metaphlan_jobs.sh
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

host_removed_sub="$HOST_REMOVED_DIR/$SUBFOLD"
output_dir="$METAPHLAN_DIR/$SUBFOLD"
log_dir="$LOG_DIR/metaphlan"
mkdir -p "$output_dir" "$log_dir"

for r1_file in "$host_removed_sub"/*_host_removed_R1.fastq.gz; do
    [ -e "$r1_file" ] || { echo "No host-removed reads found in $host_removed_sub"; exit 1; }
    base=$(basename "$r1_file" _host_removed_R1.fastq.gz)
    out_file="$output_dir/${base}_profile.txt"

    if [[ -f "$out_file" ]]; then
        echo "Skipping $base (already completed)"
        continue
    fi

    r2_file="$host_removed_sub/${base}_host_removed_R2.fastq.gz"
    comb_file="$host_removed_sub/${base}_host_removed_R1R2_combined.fastq.gz"

    sbatch --job-name="mp_${base}" \
           --cpus-per-task="$CPUS" \
           --mem="$SLURM_METAPHLAN_MEM" \
           --time="$SLURM_METAPHLAN_TIME" \
           --partition="$SLURM_PARTITION" \
           --output="$log_dir/${base}.out" \
           --error="$log_dir/${base}.err" \
           --wrap="source '$SCRIPT_DIR/common.sh'; \
                   MAMBA_EXE='$MAMBA_EXE' MAMBA_ROOT_PREFIX='$MAMBA_ROOT_PREFIX' activate_env '$METAPHLAN_ENV'; \
                   gunzip -c '$r1_file' '$r2_file' | gzip > '$comb_file'; \
                   metaphlan '$comb_file' --input_type fastq --bowtie2db '$METAPHLAN_BOWTIE2_DB' --nproc $CPUS -t rel_ab_w_read_stats -o '$out_file'; \
                   rm -f '$comb_file'"
done
