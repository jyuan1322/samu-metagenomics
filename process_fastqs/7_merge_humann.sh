#!/bin/bash
# =============================================================================
# 7_merge_humann.sh — collect all HUMAnN samples that have finished (have a
# real *_genefamilies.tsv output) and merge them into combined genefamilies,
# pathabundance, and pathcoverage tables. Safe to run repeatedly while the
# HUMAnN array job is still in progress — only picks up completed samples,
# never partial/in-progress ones.
#
# Uses a FLAT staging directory (not per-sample subdirectories) — confirmed
# that humann_join_tables does not recurse into subdirectories under --input;
# all files must sit directly in the given folder. Uses hardlinks (not
# symlinks) to avoid duplicating large files on disk; each filename already
# starts with its sample ID, so flattening doesn't cause collisions.
#
# Config: BASE_DIR, SUBFOLD, HUMANN_DIR  (see config.sh)
# Usage:  ./7_merge_humann.sh
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

output_dir="$HUMANN_DIR/$SUBFOLD"
merged_dir="$HUMANN_DIR/${SUBFOLD}_merged"
mkdir -p "$merged_dir"

staging_dir="$merged_dir/_staging"
rm -rf "$staging_dir"
mkdir -p "$staging_dir"

n_completed=0
n_incomplete=0

for sample_dir in "$output_dir"/*/; do
    [ -d "$sample_dir" ] || continue
    base=$(basename "$sample_dir")

    genefamilies_file=$(ls "$sample_dir"/*_genefamilies.tsv 2>/dev/null | head -n1 || true)
    pathabundance_file=$(ls "$sample_dir"/*_pathabundance.tsv 2>/dev/null | head -n1 || true)
    pathcoverage_file=$(ls "$sample_dir"/*_pathcoverage.tsv 2>/dev/null | head -n1 || true)

    if [[ -n "$genefamilies_file" && -n "$pathabundance_file" && -n "$pathcoverage_file" ]]; then
        ln "$(realpath "$genefamilies_file")" "$staging_dir/"
        ln "$(realpath "$pathabundance_file")" "$staging_dir/"
        ln "$(realpath "$pathcoverage_file")" "$staging_dir/"
        n_completed=$((n_completed + 1))
    else
        echo "Skipping $base (not yet complete)"
        n_incomplete=$((n_incomplete + 1))
    fi
done

echo ""
echo "Found $n_completed completed sample(s), $n_incomplete still in progress or missing output."

if [[ "$n_completed" -eq 0 ]]; then
    echo "Nothing to merge yet."
    rm -rf "$staging_dir"
    exit 0
fi

# ---------------------------------------------------------------------------
# Merge with HUMAnN's own join tool. humann_join_tables exits 0 even when it
# finds nothing to join, so each call is followed by an explicit
# non-empty-output check rather than relying on `set -e` to catch failures.
# ---------------------------------------------------------------------------
echo "Merging gene families..."
humann_join_tables --input "$staging_dir" --output "$merged_dir/joined_genefamilies.tsv" --file_name genefamilies
if [[ ! -s "$merged_dir/joined_genefamilies.tsv" ]]; then
    echo "ERROR: joined_genefamilies.tsv was not created or is empty — aborting."
    exit 1
fi

echo "Merging pathway abundance..."
humann_join_tables --input "$staging_dir" --output "$merged_dir/joined_pathabundance.tsv" --file_name pathabundance
if [[ ! -s "$merged_dir/joined_pathabundance.tsv" ]]; then
    echo "ERROR: joined_pathabundance.tsv was not created or is empty — aborting."
    exit 1
fi

echo "Merging pathway coverage..."
humann_join_tables --input "$staging_dir" --output "$merged_dir/joined_pathcoverage.tsv" --file_name pathcoverage
if [[ ! -s "$merged_dir/joined_pathcoverage.tsv" ]]; then
    echo "ERROR: joined_pathcoverage.tsv was not created or is empty — aborting."
    exit 1
fi

# ---------------------------------------------------------------------------
# Normalize gene families and pathway abundance to relative abundance, since
# raw HUMAnN output is in RPK (not comparable across samples with different
# sequencing depth).
# ---------------------------------------------------------------------------
echo "Renormalizing gene families to relative abundance..."
humann_renorm_table --input "$merged_dir/joined_genefamilies.tsv" --output "$merged_dir/joined_genefamilies_relab.tsv" --units relab

echo "Renormalizing pathway abundance to relative abundance..."
humann_renorm_table --input "$merged_dir/joined_pathabundance.tsv" --output "$merged_dir/joined_pathabundance_relab.tsv" --units relab

rm -rf "$staging_dir"

echo ""
echo "Done. Merged files in: $merged_dir"
ls -la "$merged_dir"
