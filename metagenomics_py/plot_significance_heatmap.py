"""
plot_significance_heatmap.py

Publication-quality heatmap of CLR-transformed taxon abundance for all
features significant in >=1 model (DESeq2 / lasso / FCNN / random forest),
with samples split into two column blocks (Sarc vs. NonSarc) and clustered
SEPARATELY within each block (columns are never reordered across the
Sarc/NonSarc boundary). Rows are clustered once, globally, so the same
feature order applies to both blocks. Features significant in >=2 models
(the ones that would appear in the significance grid) are marked with an
asterisk after their row label.

Significance rules (same as plot_significance_grid.py):
  DESeq2   : deseq2_padj < 0.05
  Lasso    : exp(linreg_lower_perc12_5), exp(linreg_upper_perc87_5) -> odds
             range does NOT include 1
  FCNN     : [fcnn_lower_perc12_5, fcnn_upper_perc87_5] does NOT include 0
  RF       : rf_selection_frequency > 0

Color scale: each feature (row) is z-scored across ALL included samples
(both groups combined), so color is comparable between the Sarc and
NonSarc blocks.

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
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list

from loaders import load_metagenomics

GROUP_COLORS = {"NonSarc": "#4C72B0", "Sarc": "#C44E52"}


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
    """matrix: observations x dimensions. Returns leaf order (row indices)."""
    if matrix.shape[0] < 2:
        return list(range(matrix.shape[0])), None
    Z = linkage(matrix, method=method, metric=metric)
    return leaves_list(Z), Z


def draw_dendrogram(ax, Z, n_leaves, orientation="top"):
    """Draw a dendrogram whose leaf x-positions line up with imshow column
    centers (i + 0.5), by rescaling scipy's default 10-unit leaf spacing."""
    if Z is None:
        ax.axis("off")
        return
    dd = dendrogram(Z, no_plot=True)
    for xs, ys in zip(dd["icoord"], dd["dcoord"]):
        xs_scaled = [x / 10.0 for x in xs]
        ax.plot(xs_scaled, ys, color="black", linewidth=0.7)
    ax.set_xlim(0, n_leaves)
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
    y = y.loc[mat.columns]

    # z-score each feature (row) across ALL included samples
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1, ddof=0), axis=0)
    z = z.fillna(0.0)

    # --- row clustering (global, across both groups) ---
    row_order, row_Z = cluster_order(z.values)
    z = z.iloc[row_order]
    row_labels = [clean_feature_name(f) + (" *" if f in star_features else "")
                 for f in z.index]

    # --- column clustering, separately within each group ---
    nonsarc_samples = y[y == 0].index
    sarc_samples = y[y == 1].index

    def ordered_block(samples):
        if len(samples) == 0:
            return samples, None
        sub = z[samples]
        order, Z = cluster_order(sub.T.values)  # samples x features
        return sub.columns[order], Z

    nonsarc_order, nonsarc_Z = ordered_block(nonsarc_samples)
    sarc_order, sarc_Z = ordered_block(sarc_samples)

    n_nonsarc, n_sarc = len(nonsarc_order), len(sarc_order)
    z_ordered = pd.concat([z[nonsarc_order], z[sarc_order]], axis=1)

    # --- figure layout: dendrograms on top (per block), group color bar,
    # heatmap below, blocks separated by a gap column ---
    n_feat = len(z_ordered)
    gap = max(1, round(0.02 * (n_nonsarc + n_sarc)))
    cbar_w = max(2, round(0.05 * (n_nonsarc + n_sarc)))
    fig_w = max(9, 0.18 * (n_nonsarc + n_sarc + gap + cbar_w) + 4)
    fig_h = max(6, 0.34 * n_feat + 2.5)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(3, 4, figure=fig,
                 width_ratios=[n_nonsarc, gap, n_sarc, cbar_w],
                 height_ratios=[1.1, 0.18, max(3, 0.34 * n_feat)],
                 hspace=0.03, wspace=0.02)

    ax_dendro_ns = fig.add_subplot(gs[0, 0])
    ax_dendro_s = fig.add_subplot(gs[0, 2])
    ax_group_ns = fig.add_subplot(gs[1, 0])
    ax_group_s = fig.add_subplot(gs[1, 2])
    ax_heat = fig.add_subplot(gs[2, 0:3])
    cbar_ax = fig.add_subplot(gs[2, 3])

    draw_dendrogram(ax_dendro_ns, nonsarc_Z, n_nonsarc)
    draw_dendrogram(ax_dendro_s, sarc_Z, n_sarc)

    ax_group_ns.imshow([[0] * max(n_nonsarc, 1)],
                       cmap=plt.matplotlib.colors.ListedColormap([GROUP_COLORS["NonSarc"]]),
                       aspect="auto")
    ax_group_ns.set_xticks([]); ax_group_ns.set_yticks([])
    ax_group_s.imshow([[0] * max(n_sarc, 1)],
                      cmap=plt.matplotlib.colors.ListedColormap([GROUP_COLORS["Sarc"]]),
                      aspect="auto")
    ax_group_s.set_xticks([]); ax_group_s.set_yticks([])

    vmax = np.nanpercentile(np.abs(z_ordered.values), 95) or 1.0
    im = ax_heat.imshow(z_ordered.values, aspect="auto", cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax,
                        extent=(0, n_nonsarc + gap + n_sarc, n_feat, 0))
    # blank out the gap column so it reads as a visual break between blocks
    ax_heat.axvspan(n_nonsarc, n_nonsarc + gap, color="white", zorder=3)

    ax_heat.set_yticks(np.arange(n_feat) + 0.5)
    ax_heat.set_yticklabels(row_labels, fontsize=10, fontstyle="italic")
    ax_heat.set_xticks([])
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    ax_dendro_ns.set_title("NonSarc", fontsize=13, fontweight="bold",
                           color=GROUP_COLORS["NonSarc"])
    ax_dendro_s.set_title("Sarc", fontsize=13, fontweight="bold",
                          color=GROUP_COLORS["Sarc"])

    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Row z-score (CLR abundance)", fontsize=11)
    cbar.ax.tick_params(labelsize=9)

    fig.suptitle(
        f"Features significant in \u2265{args.min_sig} model"
        f"{'s' if args.min_sig != 1 else ''} "
        f"(* = significant in \u2265{args.asterisk_min_sig} models)",
        fontsize=15, fontweight="bold", y=0.995)

    fig.savefig(f"{args.out}.pdf", bbox_inches="tight")
    fig.savefig(f"{args.out}.png", dpi=300, bbox_inches="tight")
    print(f"Saved {args.out}.pdf and {args.out}.png "
          f"({n_feat} features x {n_nonsarc + n_sarc} samples; "
          f"{len(star_features)} starred)")


if __name__ == "__main__":
    main()