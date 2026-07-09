#!/usr/bin/env python
"""
9_summary_stats_genefamilies.py — per-sample QC summary from
joined_genefamilies_relab.tsv: gene family richness (count detected) and
percent UNMAPPED, plotted as two bar charts sorted for easy outlier spotting.

Usage:
    python 9_summary_stats_genefamilies.py \
        --input /path/to/joined_genefamilies_relab.tsv \
        --output /path/to/genefamilies_summary_stats.pdf
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="joined_genefamilies_relab.tsv")
    parser.add_argument("--output", required=True, help="output image path (.pdf or .png)")
    parser.add_argument("--sample-label-len", type=int, default=20)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", index_col=0)

    # Unstratified rows only (no "|" in the ID) — total per gene family,
    # not the per-species breakdown
    unstrat = df[~df.index.str.contains(r"\|", regex=True)]

    unmapped_row = unstrat.loc[[r for r in unstrat.index if r.upper() == "UNMAPPED"]]
    features = unstrat.drop(index=unmapped_row.index, errors="ignore")

    # % UNMAPPED per sample (relab values sum to ~1 per sample, so this is
    # already a fraction — multiply by 100 for a percentage)
    pct_unmapped = (unmapped_row.sum(axis=0) * 100).sort_values(ascending=False)

    # Richness: count of gene families with non-zero relative abundance
    richness = (features > 0).sum(axis=0).sort_values(ascending=False)

    # Truncate sample labels for display
    def trunc(s, n):
        return s if len(s) <= n else s[: n - 3] + "..."

    fig, axes = plt.subplots(2, 1, figsize=(max(10, len(df.columns) * 0.25), 10))

    axes[0].bar(range(len(richness)), richness.values, color="steelblue")
    axes[0].set_xticks(range(len(richness)))
    axes[0].set_xticklabels([trunc(s, args.sample_label_len) for s in richness.index],
                             rotation=90, fontsize=7)
    axes[0].set_ylabel("Gene families detected")
    axes[0].set_title("Gene family richness per sample")

    axes[1].bar(range(len(pct_unmapped)), pct_unmapped.values, color="indianred")
    axes[1].set_xticks(range(len(pct_unmapped)))
    axes[1].set_xticklabels([trunc(s, args.sample_label_len) for s in pct_unmapped.index],
                             rotation=90, fontsize=7)
    axes[1].set_ylabel("% UNMAPPED")
    axes[1].set_title("Percent unmapped reads per sample")

    plt.tight_layout()
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Saved summary plot to {args.output}")

    # Also print a quick numeric summary to the terminal
    print("\nRichness — min/median/max:",
          richness.min(), richness.median(), richness.max())
    print("%UNMAPPED — min/median/max:",
          round(pct_unmapped.min(), 1), round(pct_unmapped.median(), 1), round(pct_unmapped.max(), 1))

    # Flag likely outliers: >2 SD from mean in either direction
    for label, series in [("richness", richness), ("%UNMAPPED", pct_unmapped)]:
        mean, std = series.mean(), series.std()
        outliers = series[(series - mean).abs() > 2 * std]
        if len(outliers) > 0:
            print(f"\nPossible {label} outliers (>2 SD from mean):")
            print(outliers)

if __name__ == "__main__":
    main()
