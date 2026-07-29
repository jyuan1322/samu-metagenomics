"""
plot_significance_heatmap.py

Publication-quality heatmap of CLR-transformed taxon abundance for all
features significant in >=1 model (DESeq2 / lasso / FCNN / random forest),
with samples split into two column blocks (NoSarc vs. Sarc) and clustered
SEPARATELY within each block (columns are never reordered across the
NoSarc/Sarc boundary). Rows are clustered once, globally, so the same
feature order applies to both blocks. Features significant in >=2 models
(the ones that would appear in the significance grid) are marked with an
asterisk after their row label.

Layout mirrors plot_nmr_quorum_results.py's Panel B heatmap exactly:
[row dendrogram | group strip + heatmap | colorbar]. Only rows get a
dendrogram; columns are still clustered within each group (NoSarc / Sarc)
but that clustering is only used to order the columns -- no column
dendrogram is drawn. The group split is shown as a single continuous
color strip with "{group} (n=...)" labels, same as the NMR script.

Significance rules (same as plot_significance_grid.py):
  DESeq2   : deseq2_padj < 0.05
  Lasso    : exp(linreg_lower_perc12_5), exp(linreg_upper_perc87_5) -> odds
             range does NOT include 1
  FCNN     : [fcnn_lower_perc12_5, fcnn_upper_perc87_5] does NOT include 0
  RF       : rf_selection_frequency > 0

Color scale: each feature (row) is z-scored across ALL included samples
(both groups combined), so color is comparable between the NoSarc and
Sarc blocks.

Requires loaders.py (from metagenomics_py) on the Python path to rebuild the
CLR-transformed abundance matrix + sarc/non-sarc labels -- no retraining
involved, just re-reading meta_filtered/meta_df and re-doing the CLR
transform.

Usage:
  python plot_significance_heatmap.py <merged_csv> \
      --input-dir /data/local/jy1008/SaMu/results/latest/metagenomics_R \
      --meta-filtered-csv meta_filtered_06242026.csv \
      --meta-df-csv meta_df_FullSaMu_06242026.csv \
      --min-sig 1 --asterisk-min-sig 2 \
      --out significance_heatmap
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list

from loaders import load_metagenomics

GROUP_COLORS = {"NoSarc": "#4DAF4A", "Sarc": "#984EA3"}
GROUP_ORDER = ["NoSarc", "Sarc"]

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]  # falls back to DejaVu Sans if Arial isn't installed
plt.rcParams["pdf.fonttype"] = 42  # editable text in Illustrator, matching plot_nmr_quorum_results.py
plt.rcParams["ps.fonttype"] = 42


def compute_significance(df):
    """Same significance rules as plot_significance_grid.py."""
    df = df.copy()
    df["sig_deseq2"] = df["deseq2_padj"] < 0.05

    lo_or = np.exp(df["linreg_lower_perc12_5"])
    hi_or = np.exp(df["linreg_upper_perc87_5"])
    df["sig_lasso"] = ~((lo_or <= 1) & (hi_or >= 1))

    lo_fc = df["fcnn_lower_perc12_5"]
    hi_fc = df["fcnn_upper_perc87_5"]
    df["sig_fcnn"] = ~((lo_fc <= 0) & (hi_fc >= 0))

    df["sig_rf"] = df["rf_selection_frequency"] > 0

    for c in ("sig_deseq2", "sig_lasso", "sig_fcnn", "sig_rf"):
        df[c] = df[c].fillna(False)

    df["n_sig"] = df[["sig_deseq2", "sig_lasso", "sig_fcnn", "sig_rf"]].sum(axis=1)
    return df


def clean_feature_name(f):
    f = str(f)
    if "s__" in f:
        f = f.split("s__")[-1]
    return f.replace("_", " ")


def cluster_order(matrix, method="average", metric="euclidean"):
    """matrix: observations x dimensions. Returns leaf order (row indices)
    and the linkage matrix (None if fewer than 2 observations)."""
    if matrix.shape[0] < 2:
        return list(range(matrix.shape[0])), None
    Z = linkage(matrix, method=method, metric=metric)
    return leaves_list(Z), Z


def draw_row_dendrogram(ax, Z):
    """Left-oriented row dendrogram -- same call as plot_nmr_quorum_results.py's
    Panel B row dendrogram (no column dendrogram is ever drawn)."""
    if Z is not None:
        dendrogram(Z, orientation="left", ax=ax, no_labels=True,
                  color_threshold=0, above_threshold_color="#555555",
                  link_color_func=lambda k: "#555555")
    ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("merged_csv")
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--meta-filtered-csv", required=True)
    ap.add_argument("--meta-df-csv", required=True)
    ap.add_argument("--min-sig", type=int, default=2,
                     help="Minimum models a feature must be significant in "
                          "to be INCLUDED in the heatmap (default: 2, i.e. "
                          "the same features that appear in the "
                          "significance grid).")
    ap.add_argument("--asterisk-min-sig", type=int, default=2,
                     help="Minimum models a feature must be significant in "
                          "to get an asterisk (default: 2, i.e. features "
                          "that would appear in the significance grid).")
    ap.add_argument("--out", default="significance_heatmap")
    args = ap.parse_args()

    merged = pd.read_csv(args.merged_csv)
    merged = compute_significance(merged)
    incl = merged[merged["n_sig"] >= args.min_sig].copy()
    if incl.empty:
        raise ValueError(f"No features significant in >={args.min_sig} model(s).")
    star_features = set(incl.loc[incl["n_sig"] >= args.asterisk_min_sig, "feature"])

    input_dir = Path(args.input_dir)
    meta_filtered = pd.read_csv(input_dir / args.meta_filtered_csv)
    meta_df = pd.read_csv(input_dir / args.meta_df_csv)

    X_clr, y, *_ = load_metagenomics(meta_filtered, meta_df, clr_transform=True)

    feats = [f for f in incl["feature"] if f in X_clr.columns]
    missing = set(incl["feature"]) - set(feats)
    if missing:
        print(f"Warning: {len(missing)} significant feature(s) not found in "
              f"the abundance matrix, skipped: {sorted(missing)[:10]}"
              f"{' ...' if len(missing) > 10 else ''}")

    mat = X_clr[feats].T  # features x samples
    y = y.loc[mat.columns]  # 0 = NoSarc, 1 = Sarc

    # z-score each feature (row) across ALL included samples
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1, ddof=0), axis=0)
    z = z.fillna(0.0)

    # --- row clustering (global, across both groups), Ward linkage --
    # same as plot_nmr_quorum_results.py's row clustering.
    row_order, row_Z = cluster_order(z.values, method="ward")
    z = z.iloc[row_order]
    row_labels = [clean_feature_name(f) + (" *" if f in star_features else "")
                 for f in z.index]

    # --- column clustering, separately within each group. Order only --
    # no column dendrogram is drawn, matching plot_nmr_quorum_results.py. ---
    sample_groups = pd.Series(np.where(y.values == 1, "Sarc", "NoSarc"),
                              index=y.index)

    def ordered_block(samples):
        if len(samples) == 0:
            return list(samples)
        sub = z[samples]
        order, _ = cluster_order(sub.T.values, method="ward")
        return list(sub.columns[order])

    col_order = []
    boundaries = []
    for g in GROUP_ORDER:
        members = [s for s in sample_groups[sample_groups == g].index
                   if s in z.columns]
        col_order.extend(ordered_block(members))
        boundaries.append(len(col_order))
    boundaries = boundaries[:-1]  # drop final boundary (end of matrix)

    ordered = z[col_order]
    n_feat = len(ordered)
    n_samp = len(col_order)

    # --- figure layout: [row dendrogram | label spacer | strip + heatmap |
    # colorbar]. The spacer column (not an axis) reserves room for the
    # larger row-label font so it doesn't overlap the dendrogram, without
    # inflating wspace globally (which would also push the colorbar away
    # from the heatmap). ---
    fig_w = max(10, 0.16 * n_samp + 5.5)
    fig_h = max(6, 0.34 * n_feat + 2.5)
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = fig.add_gridspec(1, 4, width_ratios=[0.13, 0.10, 1.0, 0.035],
                             wspace=0.05)
    dend_outer = outer[0].subgridspec(2, 1, height_ratios=[0.08, 1.0], hspace=0.02)
    heat_outer = outer[2].subgridspec(2, 1, height_ratios=[0.08, 1.0], hspace=0.02)

    ax_dend = fig.add_subplot(dend_outer[1])
    ax_strip = fig.add_subplot(heat_outer[0])
    ax_heat = fig.add_subplot(heat_outer[1])
    ax_cbar = fig.add_subplot(outer[3])

    draw_row_dendrogram(ax_dend, row_Z)

    # --- group annotation strip (top), same style as
    # plot_nmr_quorum_results.py's Panel B ---
    ax_strip.set_xlim(0, n_samp)
    ax_strip.set_ylim(0, 1)
    ax_strip.set_xticks([]); ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)
    start = 0
    for g in GROUP_ORDER:
        n = int((sample_groups.loc[col_order] == g).sum())
        ax_strip.axvspan(start, start + n, color=GROUP_COLORS[g])
        ax_strip.text(start + n / 2, 0.5, f"{g} (n={n})", ha="center", va="center",
                     fontsize=13, fontweight="bold", color="white")
        start += n

    # --- heatmap (pcolormesh, same as plot_nmr_quorum_results.py) ---
    vlim = min(3.0, np.nanmax(np.abs(ordered.values))) if ordered.size else 1.0
    im = ax_heat.pcolormesh(ordered.values, cmap="RdBu_r", vmin=-vlim, vmax=vlim,
                            edgecolors="white", linewidth=0.4)
    ax_heat.invert_yaxis()
    for b in boundaries:
        ax_heat.axvline(b, color="black", linewidth=2.2)
    ax_heat.set_xticks([])
    ax_heat.set_yticks(np.arange(n_feat) + 0.5)
    ax_heat.set_yticklabels(row_labels, fontsize=13, fontstyle="italic")
    ax_heat.tick_params(axis="y", pad=6)  # small buffer so labels don't crowd the dendrogram
    ax_heat.set_xlabel(f"Samples (n={n_samp}, clustered within group)", fontsize=12)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Row z-score (CLR abundance)", fontsize=13)
    cbar.ax.tick_params(labelsize=11)

    fig.suptitle(
        f"Features significant in \u2265{args.min_sig} model"
        f"{'s' if args.min_sig != 1 else ''} "
        f"(* = significant in \u2265{args.asterisk_min_sig} models)",
        fontsize=18, fontweight="bold", y=0.995)

    fig.savefig(f"{args.out}.pdf", bbox_inches="tight")
    fig.savefig(f"{args.out}.png", dpi=300, bbox_inches="tight")
    print(f"Saved {args.out}.pdf and {args.out}.png "
          f"({n_feat} features x {n_samp} samples; "
          f"{len(star_features)} starred)")


if __name__ == "__main__":
    main()