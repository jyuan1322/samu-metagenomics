#!/bin/bash
# =============================================================================
# 6_submit_humann_jobs.sh — submit one HUMAnN SLURM job per sample under
# $HOST_REMOVED_DIR/$SUBFOLD, using each sample's v31 MetaPhlAn profile via
# --taxonomic-profile. Requires 5b_submit_metaphlan_v31_jobs.sh to have
# completed first. SLURM version of 6_humann.sh.
#
# Config: BASE_DIR, SUBFOLD, HUMANN_*, METAPHLAN_V31_DIR, CPUS,
#         SLURM_HUMANN_*  (see config.sh)
# Usage:  ./6_submit_humann_jobs.sh
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

host_removed_sub="$HOST_REMOVED_DIR/$SUBFOLD"
profile_sub="$METAPHLAN_V31_DIR/$SUBFOLD"
output_dir="$HUMANN_DIR/$SUBFOLD"
log_dir="$LOG_DIR/humann"
mkdir -p "$output_dir" "$log_dir"

for r1_file in "$host_removed_sub"/*_host_removed_R1.fastq.gz; do
    [ -e "$r1_file" ] || { echo "No host-removed reads found in $host_removed_sub"; exit 1; }
    base=$(basename "$r1_file" _host_removed_R1.fastq.gz)
    sample_out_dir="$output_dir/$base"
    profile_file="$profile_sub/${base}_profile_v31.txt"

    if [[ -d "$sample_out_dir" ]] && ls "$sample_out_dir"/*_genefamilies.tsv >/dev/null 2>&1; then
        echo "Skipping $base (already completed)"
        continue
    fi

    if [[ ! -f "$profile_file" ]]; then
        echo "No v31 MetaPhlAn profile found for $base at $profile_file — run 5b_submit_metaphlan_v31_jobs.sh first. Skipping."
        continue
    fi

    mkdir -p "$sample_out_dir"

    r2_file="$host_removed_sub/${base}_host_removed_R2.fastq.gz"
    comb_file="$host_removed_sub/${base}_host_removed_R1R2_combined_humann.fastq.gz"

    sbatch --job-name="hu_${base}" \
           --cpus-per-task="$CPUS" \
           --mem="$SLURM_HUMANN_MEM" \
           --time="$SLURM_HUMANN_TIME" \
           --partition="$SLURM_PARTITION" \
           --output="$log_dir/${base}.out" \
           --error="$log_dir/${base}.err" \
           --wrap="source '$SCRIPT_DIR/common.sh'; \
                   MAMBA_EXE='$MAMBA_EXE' MAMBA_ROOT_PREFIX='$MAMBA_ROOT_PREFIX' activate_env '$HUMANN_ENV'; \
                   gunzip -c '$r1_file' '$r2_file' | gzip > '$comb_file'; \
                   humann --input '$comb_file' --output '$sample_out_dir' \
                          --nucleotide-database '$HUMANN_CHOCOPHLAN_DB' \
                          --protein-database '$HUMANN_UNIREF_DB' \
                          --metaphlan-options '--bowtie2db /data/bwh-comppath-seq/databases/bowtie2 --index mpa_vJun23_CHOCOPhlAnSGB_202403' \
                          --threads $CPUS; \
                   rm -f '$comb_file'"
done
