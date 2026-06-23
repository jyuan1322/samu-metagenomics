#!/bin/bash


# List of tar files (edit this list as needed)
TAR_FILES=(
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block1_anton_07022025/20231011_AV224503_4520_1-RawData-4520.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block1_anton_07022025/20231027_AV224503_4520_2-RawData-4520.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block2_anton_07212025/20240327_AV224503_4734_1-RawData-4734-attempt2.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block3_anton_07022025/20240410_AV224503_4734_2-RawData-4734.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block4_anton_07142025/0240410_AV224503_4734_3-RawData-4734.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/ugent_block4_anton_07142025/20240502AV224503_4761_1-RawData-4734.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/nextcloud_anton_06252025/20250124_AV242402_4915_1-RawData-4915.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/nextcloud_anton_06252025/20250131_AV242402_4915_2-RawData-4915.tar"
  "/data/bwh-comppath-seq/jy1008/SaMu/data/nextcloud_anton_06252025/20250131_AV242402_4915_3-RawData-4915.tar"
)

# Target directory for fastq.gz files
DATA_DIR="/data/bwh-comppath-seq/jy1008/SaMu/data/metagenomics/raw"

# Create the data directory if it doesn't exist
mkdir -p "$DATA_DIR"

# Loop through each tar file
for TAR_FILE in "${TAR_FILES[@]}"; do
  echo "Processing $TAR_FILE..."

  # Check if file exists
  if [[ ! -f "$TAR_FILE" ]]; then
    echo "  [Warning] File not found: $TAR_FILE"
    continue
  fi

  # Extract the tar file
  tar -xvf "$TAR_FILE"

  # Find and move any .fastq.gz files to the data directory
  find . -maxdepth 2 -type f -name "*.fastq.gz" -exec mv -v {} "$DATA_DIR/" \;

done

echo "Done. All .fastq.gz files have been moved to $DATA_DIR/"
