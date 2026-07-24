#!/usr/bin/env python3
"""
plot_nmr_quorum_results.py

Publication-quality two-panel figure from the outputs of 01a_nmr.R:

  Panel A: one boxplot per metabolite (individual samples jittered on top),
           split by sarcopenia status, annotated with Mann-Whitney U
           p-values (BH-adjusted across metabolites).
  Panel B: clustered heatmap of all metabolites (z-scored log10
           concentrations) — metabolites (rows) hierarchically clustered,
           samples (columns) split into groups and clustered *within* each
           group, so the group split is always visually respected. Same
           layout as the DEP2 Panel C heatmap (row dendrogram + group
           annotation strip + colorbar).

Inputs (written by 01a_nmr.R into cfg$output_dir):
  samu_nmr_FullSaMu.csv        merged raw NMR concentrations + metadata
  samu_nmr_logrel_scaled.csv   z-scored log10 concentrations (record_id x metabolite)

Usage:
  python plot_nmr_quorum_results.py \
      --raw samu_nmr_FullSaMu.csv \
      --scaled samu_nmr_logrel_scaled.csv \
      --out plot_nmr_quorum_results

  # Writes plot_nmr_quorum_results.pdf and plot_nmr_quorum_results.png (300 dpi).

  # Optional overrides:
  --group-col sarc_status_bin      # grouping column in --raw
  --group-levels NoSarc Sarc       # reference level first
  --colors "#4DAF4A" "#984EA3"     # Panel A fill colors; defaults to match
                                    # Panel B's GROUP_COLOR_CYCLE strip colors
  --top-n 10                       # Panel A: only the top 10 by p-value
                                    # (Panel B / stats CSV still show everything)
  --pseudocount 1e-6                # log10(value + X); default 1.0 matches NMR
  --value-label "Quant*Prob"        # units label (Panel A y-axis, Panel B colorbar)
  --feature-label "Species"         # row type (Panel B title); default "Metabolite"
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist

# ---------------------------------------------------------------------------
# Publication styling
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,   # editable text in Illustrator
    "ps.fonttype": 42,
    "savefig.dpi": 300,
})


# Group color cycle for the Panel B heatmap's group-annotation strip — same
# palette and same first-two-levels-in-order mapping as plot_dep2_results.py
# (NoSarc -> green, Sarc -> purple by default). Extended/repeated with a
# warning if group_levels has more than two entries.
GROUP_COLOR_CYCLE = ["#4DAF4A", "#984EA3", "#FF7F00", "#377EB8", "#E41A1C", "#A65628"]


def bh_adjust(pvals):
    """Benjamini-Hochberg FDR adjustment."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    return adj


def sig_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def load_data(raw_path, scaled_path, group_col, group_levels):
    scaled = pd.read_csv(scaled_path, index_col=0)
    metabolites = list(scaled.columns)

    raw = pd.read_csv(raw_path)
    if "record_id" not in raw.columns:
        sys.exit(f"--raw ({raw_path}) has no record_id column")
    if group_col not in raw.columns:
        sys.exit(f"--raw ({raw_path}) has no '{group_col}' column")

    missing = [m for m in metabolites if m not in raw.columns]
    if missing:
        sys.exit(
            f"{len(missing)} metabolite column(s) in --scaled are not in --raw, "
            f"e.g. {missing[:5]}. Are --raw and --scaled from the same run?"
        )

    long = raw[["record_id", group_col] + metabolites].copy()
    for m in metabolites:
        long[m] = pd.to_numeric(long[m], errors="coerce")
    long = long[long[group_col].isin(group_levels)]
    long[group_col] = pd.Categorical(long[group_col], categories=group_levels, ordered=True)

    return long, scaled.loc[scaled.index.isin(long["record_id"])], metabolites


def compute_stats(long, group_col, group_levels, metabolites):
    rows = []
    g0, g1 = group_levels
    for m in metabolites:
        a = long.loc[long[group_col] == g0, m].dropna()
        b = long.loc[long[group_col] == g1, m].dropna()
        if len(a) < 2 or len(b) < 2:
            p = np.nan
        else:
            _, p = mannwhitneyu(a, b, alternative="two-sided")
        rows.append({
            "metabolite": m,
            "n0": len(a), "n1": len(b),
            "mean0": a.mean(), "mean1": b.mean(),
            "p_value": p,
        })
    stats = pd.DataFrame(rows)
    stats["p_adj"] = bh_adjust(stats["p_value"].fillna(1.0))
    stats.loc[stats["p_value"].isna(), "p_adj"] = np.nan
    return stats.sort_values("p_value", na_position="last").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Panel A: per-metabolite boxplots
# ---------------------------------------------------------------------------
def plot_panel_a(fig, gs_cell, long, stats, group_col, group_levels, colors,
                 pseudocount=1.0, value_label="Conc."):
    """Single boxplot spanning all metabolites on the x-axis, dodged/colored
    by group — matching the original single-panel ggplot boxplot in
    01a_nmr.R (log10(value + pseudocount) transform, metabolites ordered by
    ascending raw p-value, same as R's `arrange(p_value)` reordering).
    pseudocount/value_label make this reusable for other assays with a
    different offset or units (e.g. quorum sensing uses cfg$epsilon, not
    NMR's +1, and "Quant*Prob" rather than "Conc.")."""
    ax = fig.add_subplot(gs_cell)
    rng = np.random.default_rng(0)

    metabolites = stats["metabolite"].tolist()  # already sorted by p_value
    n = len(metabolites)
    n_groups = len(group_levels)

    dodge_width = 0.8
    box_width = dodge_width / n_groups * 0.85
    offsets = np.linspace(-dodge_width / 2 + box_width / 2,
                          dodge_width / 2 - box_width / 2, n_groups)

    group_max_per_metab = []
    for i, m in enumerate(metabolites):
        group_max = -np.inf
        for gi, g in enumerate(group_levels):
            vals = np.log10(long.loc[long[group_col] == g, m].dropna().values.astype(float) + pseudocount)
            if len(vals) == 0:
                continue
            pos = i + offsets[gi]
            bp = ax.boxplot(
                [vals], positions=[pos], widths=box_width, patch_artist=True,
                showfliers=False,
                medianprops=dict(color="black", linewidth=1.0),
                boxprops=dict(linewidth=0.6), whiskerprops=dict(linewidth=0.6),
                capprops=dict(linewidth=0.6),
            )
            bp["boxes"][0].set_facecolor(colors[gi])
            bp["boxes"][0].set_alpha(0.85)
            bp["boxes"][0].set_edgecolor("black")
            group_max = max(group_max, np.nanmax(vals))

            jitter = rng.uniform(-box_width * 0.3, box_width * 0.3, size=len(vals))
            ax.scatter(pos + jitter, vals, s=3, color="black", alpha=0.4,
                      linewidths=0, zorder=3)
        group_max_per_metab.append(group_max if np.isfinite(group_max) else 0.0)

    ax.set_xticks(range(n))
    ax.set_xticklabels(metabolites, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel(f"log10({value_label} + {pseudocount:g})", fontsize=8)
    ax.set_xlim(-0.6, n - 0.4)

    y_max = max(group_max_per_metab) if group_max_per_metab else 1.0
    ax.set_ylim(top=y_max * 1.4 if y_max > 0 else 1.0)  # extra headroom for 2-line labels

    for i, row in stats.iterrows():
        if pd.notna(row["p_value"]):
            label = f"p = {row['p_value']:.2g}\n(padj = {row['p_adj']:.2g})"
        else:
            label = "n/a"
        ax.text(i, group_max_per_metab[i] + y_max * 0.03, label, ha="center", va="bottom",
                fontsize=5, linespacing=1.2)

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="black", alpha=0.85)
              for c in colors]
    ax.legend(handles, group_levels, loc="upper right", frameon=False, fontsize=7)

    return fig


# ---------------------------------------------------------------------------
# Panel B: clustered heatmap (row dendrogram + group strip), same layout as
# the DEP2 Panel C heatmap.
# ---------------------------------------------------------------------------
def cluster_order(mat, axis, metric="euclidean", method="ward"):
    """Hierarchical-cluster leaf order for rows (axis='rows') or columns
    (axis='cols') of mat. Falls back to the original order (no reordering,
    no linkage) if fewer than 3 items are present."""
    data = mat.values if axis == "rows" else mat.values.T
    labels = list(mat.index) if axis == "rows" else list(mat.columns)
    if len(labels) < 3:
        return labels, None
    dist = pdist(data, metric=metric)
    if not np.all(np.isfinite(dist)):
        dist = pdist(data, metric="euclidean")
    Z = linkage(dist, method=method)
    order = leaves_list(Z)
    return [labels[i] for i in order], Z


def cluster_columns_within_groups(mat, groups, group_order):
    """Cluster samples separately within each group level (so the group
    split is always respected), then concatenate. Returns the full ordered
    sample list and the group-boundary indices (for divider lines).

    `groups` may include samples absent from `mat` (e.g. samples with no
    NMR data at all, dropped by build_matrix() in 01a_nmr.R before the
    scaled matrix was written) — those are silently skipped here rather
    than raising, since plot_panel_b() already restricts `groups` to
    mat's columns and warns about anything dropped."""
    ordered_samples = []
    boundaries = []
    for g in group_order:
        members = [m for m in groups[groups == g].index.tolist() if m in mat.columns]
        sub = mat[members]
        order, _ = cluster_order(sub, axis="cols")
        ordered_samples.extend(order)
        boundaries.append(len(ordered_samples))
    return ordered_samples, boundaries[:-1]  # drop final boundary (end of matrix)


def plot_panel_b(fig, gs_cell, scaled, long, group_col, group_levels,
                 feature_label="Metabolite", value_label="Conc."):
    """Clustered heatmap of all features, samples grouped by status."""
    mat = scaled.T.fillna(0)  # feature x sample
    group_lookup = long.set_index("record_id")[group_col]

    # Some samples in --raw have no NMR data at all and were dropped by
    # build_matrix() in 01a_nmr.R before the scaled matrix was written, so
    # they're absent from `scaled`/`mat`. Restrict the group lookup to what
    # actually made it into the heatmap, and say so.
    n_before = len(group_lookup)
    group_lookup = group_lookup.loc[group_lookup.index.isin(mat.columns)]
    n_dropped = n_before - len(group_lookup)
    if n_dropped > 0:
        print(f"NOTE: {n_dropped} sample(s) in --raw have no NMR data in --scaled "
              f"and are excluded from the Panel B heatmap.")

    row_order, row_linkage_Z = cluster_order(mat, axis="rows", method="ward")
    col_order, col_boundaries = cluster_columns_within_groups(mat, group_lookup, group_levels)
    ordered = mat.loc[row_order, col_order]

    group_colors = {g: GROUP_COLOR_CYCLE[i % len(GROUP_COLOR_CYCLE)]
                    for i, g in enumerate(group_levels)}
    if len(group_levels) > len(GROUP_COLOR_CYCLE):
        print(f"WARNING: {len(group_levels)} group levels but only "
              f"{len(GROUP_COLOR_CYCLE)} default colors defined — colors will repeat.")

    # Layout: [row dendrogram | heatmap+strip | colorbar], matching the
    # DEP2 Panel C layout.
    outer = gs_cell.subgridspec(1, 3, width_ratios=[0.14, 1.0, 0.035], wspace=0.03)
    dend_outer = outer[0].subgridspec(2, 1, height_ratios=[0.08, 1.0], hspace=0.02)
    heat_outer = outer[1].subgridspec(2, 1, height_ratios=[0.08, 1.0], hspace=0.02)

    ax_dend = fig.add_subplot(dend_outer[1])
    ax_strip = fig.add_subplot(heat_outer[0])
    ax_heat = fig.add_subplot(heat_outer[1])
    ax_cbar = fig.add_subplot(outer[2])

    # --- row dendrogram (left) ---
    if row_linkage_Z is not None:
        dendrogram(row_linkage_Z, orientation="left", ax=ax_dend, no_labels=True,
                  color_threshold=0, above_threshold_color="#555555",
                  link_color_func=lambda k: "#555555")
    ax_dend.axis("off")

    # --- group annotation strip (top) ---
    ax_strip.set_xlim(0, len(col_order))
    ax_strip.set_ylim(0, 1)
    ax_strip.set_xticks([]); ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)
    start = 0
    for g in group_levels:
        n = int((group_lookup.loc[col_order] == g).sum())
        ax_strip.axvspan(start, start + n, color=group_colors[g])
        ax_strip.text(start + n / 2, 0.5, f"{g} (n={n})", ha="center", va="center",
                      fontsize=6, fontweight="bold", color="white")
        start += n

    # --- heatmap (pcolormesh, crisp at any zoom in the saved PDF) ---
    vlim = min(3.0, np.nanmax(np.abs(ordered.values))) if ordered.size else 1.0
    im = ax_heat.pcolormesh(ordered.values, cmap="RdBu_r", vmin=-vlim, vmax=vlim,
                            edgecolors="white", linewidth=0.4)
    ax_heat.invert_yaxis()
    for b in col_boundaries:
        ax_heat.axvline(b, color="black", linewidth=2.2)
    ax_heat.set_xticks([])
    ax_heat.set_yticks(np.arange(len(row_order)) + 0.5)
    ax_heat.set_yticklabels(row_order, fontsize=6)
    ax_heat.set_xlabel(f"Samples (n={len(col_order)}, clustered within group)", fontsize=7)

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label(f"z-score (log10 {value_label})", fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    ax_heat.set_title(f"{feature_label} levels across samples\n"
                      "(rows clustered; columns split by group, clustered within group)",
                      loc="left", fontweight="bold", fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", default="samu_nmr_FullSaMu.csv")
    ap.add_argument("--scaled", default="samu_nmr_logrel_scaled.csv")
    ap.add_argument("--group-col", default="sarc_status_bin")
    ap.add_argument("--group-levels", nargs=2, default=["NoSarc", "Sarc"])
    ap.add_argument("--colors", nargs=2, default=["#4DAF4A", "#984EA3"],
                     help="Panel A boxplot fill colors (level0, level1). Defaults to "
                          "the same GROUP_COLOR_CYCLE colors used for Panel B's "
                          "heatmap group strip.")
    ap.add_argument("--out", default="plot_nmr_quorum_results",
                     help="Output file prefix; writes <out>.pdf and <out>.png (300 dpi)")
    ap.add_argument("--stats-out", default="nmr_mannwhitney_stats.csv")
    ap.add_argument("--top-n", type=int, default=None,
                     help="Only show the top N metabolites (by ascending raw p-value) "
                          "in Panel A. Default: show all. Panel B's heatmap and "
                          "--stats-out are unaffected — always all metabolites.")
    ap.add_argument("--pseudocount", type=float, default=1.0,
                     help="Offset added before log10() in Panel A (default 1.0, "
                          "matching NMR's log10(value+1)). Quorum sensing's R script "
                          "uses cfg$epsilon (typically 1e-6) for both its boxplot and "
                          "heatmap — pass that here to match.")
    ap.add_argument("--value-label", default="Conc.",
                     help="Units label used in Panel A's y-axis and Panel B's "
                          "colorbar (default 'Conc.'; e.g. 'Quant*Prob' for quorum).")
    ap.add_argument("--feature-label", default="Metabolite",
                     help="Row/feature type label used in Panel B's title "
                          "(default 'Metabolite'; e.g. 'QSP', 'Species', "
                          "'Microbial target' for quorum).")
    args = ap.parse_args()

    long, scaled, metabolites = load_data(
        args.raw, args.scaled, args.group_col, args.group_levels
    )
    stats = compute_stats(long, args.group_col, args.group_levels, metabolites)
    stats.to_csv(args.stats_out, index=False)

    stats_panel_a = stats.head(args.top_n) if args.top_n else stats

    n = len(stats_panel_a)
    panel_a_h = 5.0
    panel_b_h = 7.5
    fig_h = panel_a_h + panel_b_h
    fig_w = max(11, n * 0.9)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(2, 1, height_ratios=[panel_a_h, panel_b_h], hspace=0.5)

    plot_panel_a(fig, gs[0], long, stats_panel_a, args.group_col, args.group_levels, args.colors,
                pseudocount=args.pseudocount, value_label=args.value_label)
    plot_panel_b(fig, gs[1], scaled, long, args.group_col, args.group_levels,
                feature_label=args.feature_label, value_label=args.value_label)

    # fig.text(0.01, 0.99, "A", fontsize=14, fontweight="bold", va="top")
    # fig.text(0.01, panel_b_h / fig_h - 0.02, "B", fontsize=14, fontweight="bold", va="top")

    # --out is a prefix, not a full filename — strip an accidental .pdf/.png
    # extension so users passing the old-style full filename still work.
    prefix = args.out
    if prefix.lower().endswith((".pdf", ".png")):
        prefix = prefix.rsplit(".", 1)[0]

    fig.savefig(f"{prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight")
    print(f"Wrote {prefix}.pdf")
    print(f"Wrote {prefix}.png")
    print(f"Wrote {args.stats_out}")


if __name__ == "__main__":
    main()