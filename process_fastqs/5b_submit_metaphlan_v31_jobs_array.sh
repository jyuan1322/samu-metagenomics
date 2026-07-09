#!/bin/bash
# =============================================================================
# 5b_submit_metaphlan_v31_jobs.sh — submit ONE SLURM array job covering every
# sample under $HOST_REMOVED_DIR/$SUBFOLD, running MetaPhlAn (v31, pre-SGB)
# against each. Replaces the old per-sample sbatch-in-a-loop approach, which
# was prone to "Slurm temporarily unable to accept job" failures under load
# from firing many sbatch calls in rapid succession — a single array
# submission avoids that entirely.
#
# Each array task looks up its own sample from a generated sample list via
# $SLURM_ARRAY_TASK_ID, so no per-sample --wrap string is built at submission
# time.
#
# Config: BASE_DIR, SUBFOLD, METAPHLAN_V31_*, CPUS, SLURM_METAPHLAN_V31_*
#         (see config.sh)
# Usage:  ./5b_submit_metaphlan_v31_jobs.sh
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

host_removed_sub="$HOST_REMOVED_DIR/$SUBFOLD"
output_dir="$METAPHLAN_V31_DIR/$SUBFOLD"
log_dir="$LOG_DIR/metaphlan_v31"
mkdir -p "$output_dir" "$log_dir"

# Cap on concurrently-running array tasks, to avoid hogging the partition.
# Override at runtime with: ARRAY_THROTTLE=5 ./5b_submit_metaphlan_v31_jobs.sh
ARRAY_THROTTLE="${ARRAY_THROTTLE:-10}"

# ---------------------------------------------------------------------------
# Build the sample list — only samples that don't already have output.
# ---------------------------------------------------------------------------
sample_list="$log_dir/sample_list.txt"
: > "$sample_list"

for r1_file in "$host_removed_sub"/*_host_removed_R1.fastq.gz; do
    [ -e "$r1_file" ] || { echo "No host-removed reads found in $host_removed_sub"; exit 1; }
    base=$(basename "$r1_file" _host_removed_R1.fastq.gz)
    out_file="$output_dir/${base}_profile_v31.txt"

    if [[ -f "$out_file" ]]; then
        echo "Skipping $base (already completed)"
        continue
    fi
    echo "$base" >> "$sample_list"
done

n_samples=$(wc -l < "$sample_list")

if [[ "$n_samples" -eq 0 ]]; then
    echo "Nothing to do — all samples already have output."
    exit 0
fi

echo "Submitting array job for $n_samples sample(s) (throttled to $ARRAY_THROTTLE concurrent tasks)"

# ---------------------------------------------------------------------------
# Submit the array. Each task reads its own sample from $sample_list via
# $SLURM_ARRAY_TASK_ID (1-indexed, matching `sed -n Np`).
# ---------------------------------------------------------------------------
sbatch --job-name="mpv31_array" \
       --array=1-"${n_samples}%${ARRAY_THROTTLE}" \
       --cpus-per-task="$CPUS" \
       --mem="$SLURM_METAPHLAN_V31_MEM" \
       --time="$SLURM_METAPHLAN_V31_TIME" \
       --partition="$SLURM_PARTITION" \
       --output="$log_dir/%a.out" \
       --error="$log_dir/%a.err" \
       --wrap="source '$SCRIPT_DIR/common.sh'; \
               base=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" '$sample_list'); \
               echo \"Array task \$SLURM_ARRAY_TASK_ID -> sample \$base\"; \
               MAMBA_EXE='$MAMBA_EXE' MAMBA_ROOT_PREFIX='$MAMBA_ROOT_PREFIX' activate_env '$METAPHLAN_V31_ENV'; \
               r1_file='$host_removed_sub'/\${base}_host_removed_R1.fastq.gz; \
               r2_file='$host_removed_sub'/\${base}_host_removed_R2.fastq.gz; \
               comb_file='$host_removed_sub'/\${base}_host_removed_R1R2_combined_v31.fastq.gz; \
               out_file='$output_dir'/\${base}_profile_v31.txt; \
               gunzip -c \"\$r1_file\" \"\$r2_file\" | gzip > \"\$comb_file\"; \
               mp_max_retries=3; \
               mp_attempt=1; \
               until metaphlan \"\$comb_file\" --input_type fastq --bowtie2db '$METAPHLAN_V31_BOWTIE2_DB' --index '$METAPHLAN_V31_INDEX' --nproc $CPUS -t rel_ab_w_read_stats -o \"\$out_file\"; do \
                   rm -f \"\$out_file\"; \
                   if [[ \$mp_attempt -ge \$mp_max_retries ]]; then echo 'MetaPhlAn failed after retries'; rm -f \"\$comb_file\"; exit 1; fi; \
                   echo \"MetaPhlAn failed, attempt \$mp_attempt/\$mp_max_retries, retrying...\"; \
                   sleep 30; \
                   mp_attempt=\$((mp_attempt + 1)); \
               done; \
               rm -f \"\$comb_file\""

echo ""
echo "Array submitted. Logs will land in $log_dir/<task_id>.out / .err"
echo "Check status with: squeue -u \$USER"
echo "Check for failures with: sacct -j <jobid> --format=JobID,JobName,State,ExitCode"
