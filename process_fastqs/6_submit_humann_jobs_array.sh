#!/bin/bash
# =============================================================================
# 6_submit_humann_jobs.sh — submit ONE SLURM array job running HUMAnN for
# every sample that already has a completed v31 MetaPhlAn profile
# (from 5b_submit_metaphlan_v31_jobs.sh) and doesn't yet have HUMAnN output.
# Safe to rerun repeatedly as more v31 profiles finish — it only submits for
# samples with a profile ready right now.
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

ARRAY_THROTTLE="${ARRAY_THROTTLE:-10}"

# ---------------------------------------------------------------------------
# Build the sample list from EXISTING v31 profiles only — not from the full
# raw-read sample set. Skips samples whose HUMAnN output already exists.
# ---------------------------------------------------------------------------
sample_list="$log_dir/sample_list.txt"
: > "$sample_list"

shopt -s nullglob
profile_files=("$profile_sub"/*_profile_v31.txt)
shopt -u nullglob

if [[ ${#profile_files[@]} -eq 0 ]]; then
    echo "No v31 MetaPhlAn profiles found yet in $profile_sub"
    exit 0
fi

for profile_file in "${profile_files[@]}"; do
    base=$(basename "$profile_file" _profile_v31.txt)
    sample_out_dir="$output_dir/$base"

    if [[ -d "$sample_out_dir" ]] && ls "$sample_out_dir"/*_genefamilies.tsv >/dev/null 2>&1; then
        echo "Skipping $base (HUMAnN already completed)"
        continue
    fi
    echo "$base" >> "$sample_list"
done

n_samples=$(wc -l < "$sample_list")

if [[ "$n_samples" -eq 0 ]]; then
    echo "Nothing to do — all samples with a v31 profile already have HUMAnN output."
    exit 0
fi

echo "Submitting HUMAnN array job for $n_samples sample(s) (throttled to $ARRAY_THROTTLE concurrent tasks)"

# ---------------------------------------------------------------------------
# Submit the array.
# ---------------------------------------------------------------------------
sbatch --job-name="humann_array" \
       --array=1-"${n_samples}%${ARRAY_THROTTLE}" \
       --cpus-per-task="$CPUS" \
       --mem="$SLURM_HUMANN_MEM" \
       --time="$SLURM_HUMANN_TIME" \
       --partition="$SLURM_PARTITION" \
       --output="$log_dir/%a.out" \
       --error="$log_dir/%a.err" \
       --wrap="source '$SCRIPT_DIR/common.sh'; \
               base=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" '$sample_list'); \
               echo \"Array task \$SLURM_ARRAY_TASK_ID -> sample \$base\"; \
               MAMBA_EXE='$MAMBA_EXE' MAMBA_ROOT_PREFIX='$MAMBA_ROOT_PREFIX' activate_env '$HUMANN_ENV'; \
               r1_file='$host_removed_sub'/\${base}_host_removed_R1.fastq.gz; \
               r2_file='$host_removed_sub'/\${base}_host_removed_R2.fastq.gz; \
               comb_file='$host_removed_sub'/\${base}_host_removed_R1R2_combined_humann.fastq.gz; \
               profile_file='$profile_sub'/\${base}_profile_v31.txt; \
               sample_out_dir='$output_dir'/\${base}; \
               mkdir -p \"\$sample_out_dir\"; \
               gunzip -c \"\$r1_file\" \"\$r2_file\" | gzip > \"\$comb_file\"; \
               humann --input \"\$comb_file\" --output \"\$sample_out_dir\" \
                      --taxonomic-profile \"\$profile_file\" \
                      --nucleotide-database '$HUMANN_CHOCOPHLAN_DB' \
                      --protein-database '$HUMANN_UNIREF_DB' \
                      --threads $CPUS \
                      --remove-temp-output; \
               rm -f \"\$comb_file\""

echo ""
echo "Array submitted. Logs will land in $log_dir/<task_id>.out / .err"
echo "Check status with: squeue -u \$USER"
