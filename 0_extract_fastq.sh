#!/bin/bash

TAR_FILE="$1"
OUTDIR="$2"

echo "Extracting: $TAR_FILE into $OUTDIR"

# Extract the tar file (preserve the original tar)
tar -xvf "$TAR_FILE" -C "$OUTDIR"

# Move all .fastq.gz files found to OUTDIR root
find "$OUTDIR" -type f -name "*.fastq.gz" -exec mv -v {} "$OUTDIR" \;

echo "Finished processing $TAR_FILE"
