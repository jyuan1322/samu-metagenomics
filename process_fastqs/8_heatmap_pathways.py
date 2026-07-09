#!/usr/bin/env python
"""
8_heatmap_pathways.py — clustered heatmap of the top N pathways by mean
relative abundance across samples, from joined_pathabundance_relab.tsv.

Usage:
    python 8_heatmap_pathways.py \
        --input /path/to/joined_pathabundance_relab.tsv \
        --output /path/to/pathway_heatmap.png \
        --top-n 30
"""
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="joined_pathabundance_relab.tsv")
    parser.add_argument("--output", required=True, help="output image path (.png or .pdf)")
    parser.add_argument("--top-n", type=int, default=30, help="number of top pathways to show")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", index_col=0)

    # Keep only unstratified rows (pathway totals, not species-level breakdown)
    # Stratified rows contain "|" in the pathway ID (e.g. "PWY123|g__Species...")
    df = df[~df.index.str.contains(r"\|", regex=True)]

    # Drop the UNMAPPED/UNINTEGRATED bookkeeping rows — not informative for a
    # pathway-composition heatmap
    drop_rows = [r for r in df.index if r.upper() in ("UNMAPPED", "UNINTEGRATED")]
    df = df.drop(index=drop_rows, errors="ignore")

    # Select top N pathways by mean relative abundance across samples
    top_pathways = df.mean(axis=1).sort_values(ascending=False).head(args.top_n).index
    df_top = df.loc[top_pathways]

    # Log-scale for visualization — relative abundance data is heavily
    # right-skewed, and raw-scale heatmaps tend to be dominated by a few
    # very abundant pathways with everything else washed out
    df_log = np.log10(df_top + 1e-6)

    # Shorten pathway labels (MetaCyc IDs are often long descriptive strings)
    # df_log.index = [i if len(i) <= 60 else i[:57] + "..." for i in df_log.index]

    # Sample (column) labels: truncate to 20 chars
    df_log.columns = [
        c if len(c) <= 20 else c[:17] + "..."
        for c in df_log.columns
    ]

    plt.figure(figsize=(max(10, df_log.shape[1] * 0.3), max(8, args.top_n * 0.25)))
    g = sns.clustermap(
        df_log,
        cmap="viridis",
        figsize=(max(10, df_log.shape[1] * 0.3), max(8, args.top_n * 0.25)),
        xticklabels=True,
        yticklabels=True,
        cbar_pos=(1.005, 0.3, 0.02, 0.4),   # move legend out of frame
        cbar_kws={"label": "log10(relative abundance)"},
    )
    g.ax_heatmap.set_xlabel("Sample")
    g.ax_heatmap.set_ylabel("Pathway")
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=6)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=7)

    g.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Saved heatmap to {args.output}")
    print(f"Included {len(df_log)} pathways x {df_log.shape[1]} samples")

if __name__ == "__main__":
    main()
