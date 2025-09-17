#!/bin/bash

# List of tar files
TAR_FILES=(
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block1_anton_07022025/20231027_AV224503_4520_2-RawData-4520.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block2_anton_07212025/20240327_AV224503_4734_1-RawData-4734.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block3_anton_07022025/20240410_AV224503_4734_2-RawData-4734.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block4_anton_07142025/0240410_AV224503_4734_3-RawData-4734.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block4_anton_07142025/20240502AV224503_4761_1-RawData-4734.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/nextcloud_anton_06252025/20250124_AV242402_4915_1-RawData-4915.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/nextcloud_anton_06252025/20250131_AV242402_4915_2-RawData-4915.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/nextcloud_anton_06252025/20250131_AV242402_4915_3-RawData-4915.tar"
)

# Output base directory
DATA_DIR="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics/raw"
mkdir -p "$DATA_DIR"

# Submit a job for each TAR file
for TAR_FILE in "${TAR_FILES[@]}"; do
  BASENAME=$(basename "$TAR_FILE" .tar)
  OUTDIR="$DATA_DIR/$BASENAME"
  mkdir -p "$OUTDIR"

  cat <<EOF > "$OUTDIR/job.sh"
#!/bin/bash
#SBATCH --job-name=untar_${BASENAME}
#SBATCH --output=$OUTDIR/slurm-%j.out
#SBATCH --error=$OUTDIR/slurm-%j.err
#SBATCH --partition=bwh_comppath
#SBATCH --time=3:00:00
#SBATCH --mem=20G

bash extract_fastq.sh "$TAR_FILE" "$OUTDIR"
EOF

  sbatch "$OUTDIR/job.sh"
done

echo "All SLURM jobs submitted."
