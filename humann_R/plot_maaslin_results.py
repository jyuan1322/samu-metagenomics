#!/usr/bin/env python3
"""
plot_maaslin_results.py

Three-panel figure summarizing a MaAsLin3 all_results.tsv for a single
metadata term of interest (e.g. sarc_status_bin):

  Panel A: volcano plot (abundance-model coefficient vs. -log10(p-value)),
           colored by which sub-model (abundance vs. prevalence) drove the
           joint significance for that feature.
  Panel B: horizontal bar chart of the top-N features by q-value, signed by
           effect direction, colored the same way as Panel A, with p- and
           q-values reported to the right of each bar.
  Panel C (optional): clustered heatmap of per-sample enrichment (row
           z-scored log-abundance) for the same top-N pathways shown in
           Panel B. Pathways (rows) are hierarchically clustered; samples
           (columns) are first split into groups by --group-col (e.g.
           sarc_status_bin), then hierarchically clustered *within* each
           group, so the group split is always visually respected and
           clustering only reorders samples inside a group.

Panel C requires two additional inputs beyond all_results.tsv, since
per-sample abundance values and group membership aren't present in
MaAsLin3's own results table: --abundance-table (the filtered feature
matrix used as MaAsLin3 input) and --sample-metadata (sample -> group
mapping). See the module-level NOTE below for how to export these from the
R pipeline. If neither is supplied, the script falls back to the original
two-panel figure.

Why q-values are recomputed here rather than trusting MaAsLin3's own
qval_joint column directly: qval_joint is FDR-corrected across every
metadata term tested (all fixed effects, not just the one of interest),
which makes it needlessly conservative for a single-variable figure. This
script isolates the metadata term with --metadata, then BH-corrects just
that term's pval_joint across features — matching the guidance in
MaAsLin3's own documentation ("it may be preferable to FDR correct just the
p-values from the variables of interest").

NOTE — exporting the two Panel C input files from R (e.g. in 03_maaslin.R,
after feature_mat / meta_aligned are loaded, before running maaslin3):

    write.csv(
      data.frame(feature = rownames(feature_mat), feature_mat, check.names = FALSE),
      tag_filename("filtered_features_abundance.csv"), row.names = FALSE
    )
    write.csv(
      meta_aligned[, c("File_ID", "sarc_status_bin")],
      tag_filename("sample_metadata_for_heatmap.csv"), row.names = FALSE
    )

Usage:
    python plot_maaslin_results.py \
        --input all_results.tsv \
        --metadata sarc_status_bin \
        --qval-threshold 0.1 \
        --top-n 15 \
        --abundance-table filtered_features_abundance.csv \
        --sample-metadata sample_metadata_for_heatmap.csv \
        --sample-id-col File_ID \
        --group-col sarc_status_bin \
        --output maaslin_volcano_top_features

Produces <output>.pdf and <output>.png (300 DPI), plus <output>_table.csv
with the recomputed q-values for the isolated metadata term, so the numbers
behind the figure are auditable rather than baked silently into a plot.
"""

import argparse
import sys
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist


# -----------------------------------------------------------------------------
# BH (Benjamini-Hochberg) correction, implemented directly rather than via
# statsmodels/scipy, so this script has no dependency beyond
# pandas/numpy/matplotlib/scipy. Matches R's p.adjust(method = "BH").
# -----------------------------------------------------------------------------
def bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


# -----------------------------------------------------------------------------
# Load + reshape a MaAsLin3 all_results.tsv down to one row per feature for
# the metadata term of interest, with:
#   - abund_coef / abund_pval  : from the "abundance" model row
#   - prev_coef  / prev_pval   : from the "prevalence" model row (may be NA
#                                  if that feature had a degenerate/error case)
#   - driver                   : "abundance" or "prevalence", whichever
#                                  sub-model had the smaller individual
#                                  p-value (falls back to whichever model
#                                  didn't error, if only one succeeded)
#   - qval_isolated             : BH-corrected pval_joint, restricted to just
#                                  this metadata term (see module docstring)
# -----------------------------------------------------------------------------
def prepare_feature_table(df: pd.DataFrame, metadata: str) -> pd.DataFrame:
    sub = df[df["metadata"] == metadata].copy()
    if sub.empty:
        raise ValueError(
            f"No rows found with metadata == '{metadata}'. "
            f"Available metadata values: {sorted(df['metadata'].unique())}"
        )

    abund = sub[sub["model"] == "abundance"].set_index("feature")
    prev = sub[sub["model"] == "prevalence"].set_index("feature")

    features = sorted(set(abund.index) | set(prev.index))
    rows = []
    for feat in features:
        a = abund.loc[feat] if feat in abund.index else None
        p = prev.loc[feat] if feat in prev.index else None

        a_ok = a is not None and pd.isna(a.get("error"))
        p_ok = p is not None and pd.isna(p.get("error"))

        if not a_ok and not p_ok:
            continue  # both sub-models failed for this feature; nothing to plot

        # pval_joint is identical across the abundance/prevalence rows for a
        # given feature (it's already the combined test) — pull it from
        # whichever row succeeded.
        pval_joint = a["pval_joint"] if a_ok else p["pval_joint"]

        a_pval_ind = a["pval_individual"] if a_ok else np.nan
        p_pval_ind = p["pval_individual"] if p_ok else np.nan

        if a_ok and p_ok:
            driver = "abundance" if a_pval_ind <= p_pval_ind else "prevalence"
        elif a_ok:
            driver = "abundance"
        else:
            driver = "prevalence"

        rows.append({
            "feature": feat,
            "abund_coef": a["coef"] if a_ok else np.nan,
            "prev_coef": p["coef"] if p_ok else np.nan,
            "pval_joint": pval_joint,
            "driver": driver,
            "N_not_zero": (a["N_not_zero"] if a_ok else p["N_not_zero"]),
        })

    out = pd.DataFrame(rows)
    out["qval_isolated"] = bh_qvalues(out["pval_joint"].values)
    return out.sort_values("pval_joint").reset_index(drop=True)


# -----------------------------------------------------------------------------
# omit_pathways — drop rows whose feature matches an entry in omit_list.
# Matches on the pathway ID portion (text before the first ":", e.g.
# "PWY-101" from "PWY-101: some pathway name") OR the full feature string,
# so either a bare ID or a full copy-pasted feature name works as an omit
# entry. Case-sensitive, exact match only (not substring) — avoids
# accidentally dropping unrelated pathways that happen to share a substring.
# -----------------------------------------------------------------------------
def omit_pathways(df: pd.DataFrame, omit_list: list) -> pd.DataFrame:
    if not omit_list:
        return df
    omit_set = set(omit_list)
    pathway_ids = df["feature"].str.split(":", n=1).str[0].str.strip()
    keep_mask = ~(df["feature"].isin(omit_set) | pathway_ids.isin(omit_set))
    n_dropped = (~keep_mask).sum()
    if n_dropped > 0:
        dropped = df.loc[~keep_mask, "feature"].tolist()
        print(f"Omitting {n_dropped} pathway(s) per --omit-pathways/--omit-file: "
              f"{dropped}")
    unmatched = omit_set - set(pathway_ids) - set(df["feature"])
    if unmatched:
        print(f"WARNING: {len(unmatched)} omit entry(ies) did not match any "
              f"feature in the input and had no effect: {sorted(unmatched)}")
    return df[keep_mask].reset_index(drop=True)


# -----------------------------------------------------------------------------
# Cosmetic helpers
# -----------------------------------------------------------------------------
def clean_feature_label(feature: str, max_len: int = 45) -> str:
    """Strip a leading 'PWY-XXXX: ' style ID prefix for display, keep the
    readable name, truncate long names for axis/label space."""
    if ":" in feature:
        label = feature.split(":", 1)[1].strip()
    else:
        label = feature
    if len(label) > max_len:
        label = label[: max_len - 1] + "\u2026"
    return label


def wrap_feature_label(feature: str, width: int = 28) -> str:
    """Like clean_feature_label, but wraps onto multiple lines (via
    textwrap.fill) instead of truncating with an ellipsis — for contexts
    (e.g. the bar chart's y-axis) where the full name is preferred and
    vertical space for extra lines is available."""
    if ":" in feature:
        label = feature.split(":", 1)[1].strip()
    else:
        label = feature
    return textwrap.fill(label, width=width)


DRIVER_COLORS = {"abundance": "#2166AC", "prevalence": "#B2182B"}
DRIVER_LABELS = {"abundance": "Abundance-driven", "prevalence": "Prevalence-driven"}

# Default group color cycle for Panel C's group-annotation strip. Extended
# as needed if --group-col has more than 4 levels (falls back to a warning).
GROUP_COLOR_CYCLE = ["#4DAF4A", "#984EA3", "#FF7F00", "#377EB8", "#E41A1C", "#A65628"]


def set_publication_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,   # embed fonts as editable text in PDF, not outlines
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# -----------------------------------------------------------------------------
# Panel A: volcano plot
# -----------------------------------------------------------------------------
def compute_bh_pvalue_threshold(feat_df: pd.DataFrame, qval_threshold: float):
    """Convert a q-value significance threshold into a single p-value line for
    a raw-p volcano plot. Two cases, since BH significance is rank-dependent
    rather than a fixed p-value cutoff:

      - If >=1 feature already clears qval_threshold: the true empirical BH
        cutoff is the LARGEST pval_joint among those significant features —
        return that, with is_empirical=True.
      - If none do: there's no empirical cutoff to report. Fall back to the
        rank-1 BH critical value (alpha / n) — the smallest, strictest
        boundary, representing "how small the single best p-value would have
        needed to be." Returned with is_empirical=False so the plot can label
        it as a conservative reference rather than an observed cutoff.
    """
    n = len(feat_df)
    sig = feat_df[feat_df["qval_isolated"] < qval_threshold]
    if len(sig) > 0:
        return sig["pval_joint"].max(), True
    return qval_threshold / n, False


def plot_volcano(ax, feat_df: pd.DataFrame, qval_threshold: float, label_top_n: int):
    x = feat_df["abund_coef"].values
    y = -np.log10(feat_df["pval_joint"].values.clip(min=1e-300))

    for driver, sub in feat_df.groupby("driver"):
        idx = feat_df.index.isin(sub.index)
        ax.scatter(
            x[idx], y[idx],
            s=22, linewidths=0.4, edgecolors="white",
            color=DRIVER_COLORS[driver], alpha=0.85,
            label=DRIVER_LABELS[driver], zorder=3,
        )

    p_threshold, is_empirical = compute_bh_pvalue_threshold(feat_df, qval_threshold)
    sig_line_y = -np.log10(p_threshold)
    ax.axhline(sig_line_y, color="grey", linestyle="--", linewidth=0.7, zorder=1)
    ax.axvline(0, color="grey", linestyle="-", linewidth=0.6, zorder=1)
    line_label = (
        f"p = {p_threshold:.2g}\n(q = {qval_threshold} cutoff)"
        if is_empirical else
        f"p = {p_threshold:.2g}\n(rank-1 threshold; none passed)"
    )
    ax.text(0.98, sig_line_y, line_label, transform=ax.get_yaxis_transform(),
            va="bottom", ha="right", fontsize=6, color="grey", linespacing=1.3)

    top = feat_df.nsmallest(label_top_n, "qval_isolated")
    for _, row in top.iterrows():
        ax.annotate(
            clean_feature_label(row["feature"], max_len=28),
            xy=(row["abund_coef"], -np.log10(max(row["pval_joint"], 1e-300))),
            xytext=(4, 3), textcoords="offset points",
            fontsize=6, color="black",
        )

    ax.set_xlabel("Abundance-model coefficient (log$_2$ fold change)")
    ax.set_ylabel(r"$-\log_{10}(p\mathrm{-value})$")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper left", handletextpad=0.4, borderaxespad=0.2)


# -----------------------------------------------------------------------------
# Panel B: top-N features, horizontal bar chart
# -----------------------------------------------------------------------------
def plot_top_features(ax, top: pd.DataFrame, top_n: int, qval_threshold: float,
                      label_wrap_width: int = 28):
    labels = [wrap_feature_label(f, width=label_wrap_width) for f in top["feature"]]
    colors = [DRIVER_COLORS[d] for d in top["driver"]]
    values = top["abund_coef"].values

    max_lines = max(label.count("\n") + 1 for label in labels)
    row_spacing = 1.0 + 0.45 * (max_lines - 1)
    ypos = np.arange(len(top)) * row_spacing
    bar_height = 0.65 * min(row_spacing, 1.6)

    ax.barh(ypos, values, color=colors, height=bar_height, edgecolor="white", linewidth=0.4)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_ylim(-row_spacing * 0.75, ypos[-1] + row_spacing * 0.75)
    ax.axvline(0, color="grey", linewidth=0.6)

    for y, (_, row) in zip(ypos, top.iterrows()):
        sig_marker = "*" if row["qval_isolated"] < qval_threshold else ""
        ax.text(1.03, y, f"p={row['pval_joint']:.2g}, q={row['qval_isolated']:.2g}{sig_marker}",
                transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=6, color="black")

    ax.set_xlabel("Abundance-model coefficient (log$_2$ fold change)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"Top {top_n} pathways by q-value", loc="left", fontweight="bold", fontsize=8)


# -----------------------------------------------------------------------------
# Panel C: clustered heatmap of per-sample enrichment for the top pathways
# -----------------------------------------------------------------------------
def load_abundance_table(path: str) -> pd.DataFrame:
    """Feature x sample table: first column 'feature' (or unnamed -> treated
    as the index), remaining columns are samples. Handles CSV or TSV by
    extension."""
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    id_col = "feature" if "feature" in df.columns else df.columns[0]
    df = df.set_index(id_col)
    df.index.name = "feature"
    return df


def load_sample_metadata(path: str, sample_id_col: str, group_col: str) -> pd.DataFrame:
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    missing = {sample_id_col, group_col} - set(df.columns)
    if missing:
        sys.exit(f"--sample-metadata is missing expected column(s): {missing}. "
                 f"Found columns: {list(df.columns)}")
    df = df.set_index(sample_id_col)
    df.index = df.index.astype(str)
    return df[[group_col]]


def build_enrichment_matrix(abundance_df: pd.DataFrame, meta_df: pd.DataFrame,
                            top_features: list, group_col: str):
    """Restrict to the top pathways x samples with known group membership,
    then row-z-score log10(abundance + pseudocount) — each pathway's values
    are expressed relative to its own mean/sd across samples ("enrichment"),
    which is what makes a cross-pathway heatmap color scale meaningful (raw
    relative abundances span many orders of magnitude across pathways and
    would otherwise all look like a single hot/cold row)."""
    abundance_df = abundance_df.copy()
    abundance_df.columns = abundance_df.columns.astype(str)

    missing_features = [f for f in top_features if f not in abundance_df.index]
    if missing_features:
        sys.exit(
            "The following top-panel-B feature(s) were not found in "
            "--abundance-table (row identifiers must match all_results.tsv's "
            f"'feature' column exactly): {missing_features}"
        )

    shared_samples = [s for s in abundance_df.columns if s in meta_df.index]
    n_dropped_abund = abundance_df.shape[1] - len(shared_samples)
    n_dropped_meta = meta_df.shape[0] - len(shared_samples)
    if n_dropped_abund > 0:
        print(f"NOTE: {n_dropped_abund} sample(s) in --abundance-table have no "
              f"matching row in --sample-metadata and are excluded from Panel C.")
    if n_dropped_meta > 0:
        print(f"NOTE: {n_dropped_meta} sample(s) in --sample-metadata have no "
              f"matching column in --abundance-table and are excluded from Panel C.")
    if len(shared_samples) < 2:
        sys.exit("Fewer than 2 samples matched between --abundance-table and "
                 "--sample-metadata — cannot build Panel C.")

    mat = abundance_df.loc[top_features, shared_samples].astype(float)
    groups = meta_df.loc[shared_samples, group_col]

    pseudocount = mat[mat > 0].min().min()
    if pd.isna(pseudocount) or pseudocount <= 0:
        pseudocount = 1e-6
    log_mat = np.log10(mat + pseudocount)

    row_mean = log_mat.mean(axis=1)
    row_std = log_mat.std(axis=1, ddof=0)
    # Guard against a (rare, post-filtering) zero-variance row: z-score would
    # be 0/0. Leave it as a flat 0 (no enrichment signal) rather than NaN,
    # so it doesn't break clustering distances.
    row_std_safe = row_std.replace(0, np.nan)
    z_mat = log_mat.sub(row_mean, axis=0).div(row_std_safe, axis=0).fillna(0.0)

    return z_mat, groups


def cluster_order(mat: pd.DataFrame, axis: str, metric: str = "euclidean",
                  method: str = "average"):
    """Hierarchical-cluster leaf order for rows (axis='rows') or columns
    (axis='cols') of mat. Falls back to the original order (no reordering,
    no linkage) if fewer than 3 items are present, since a dendrogram over
    1-2 items is degenerate/uninformative and pdist can error on n<2."""
    data = mat.values if axis == "rows" else mat.values.T
    labels = list(mat.index) if axis == "rows" else list(mat.columns)
    if len(labels) < 3:
        return labels, None
    dist = pdist(data, metric=metric)
    if not np.all(np.isfinite(dist)):
        # e.g. correlation distance undefined for a constant row/column pair
        # — fall back to euclidean, which is always finite for finite input.
        dist = pdist(data, metric="euclidean")
    Z = linkage(dist, method=method)
    order = leaves_list(Z)
    return [labels[i] for i in order], Z


def cluster_columns_within_groups(z_mat: pd.DataFrame, groups: pd.Series, group_order: list):
    """Cluster samples separately within each group level (so the group
    split is always respected), then concatenate. Returns the full ordered
    sample list and the group-boundary indices (for divider lines)."""
    ordered_samples = []
    boundaries = []
    for g in group_order:
        members = groups[groups == g].index.tolist()
        sub = z_mat[members]
        order, _ = cluster_order(sub, axis="cols")
        ordered_samples.extend(order)
        boundaries.append(len(ordered_samples))
    return ordered_samples, boundaries[:-1]  # drop final boundary (end of matrix)


def plot_heatmap(fig, gs_cell, z_mat: pd.DataFrame, groups: pd.Series,
                 group_order: list, label_wrap_width: int = 28):
    row_order, row_linkage_Z = cluster_order(z_mat, axis="rows")
    col_order, col_boundaries = cluster_columns_within_groups(z_mat, groups, group_order)
    ordered = z_mat.loc[row_order, col_order]

    group_colors = {g: GROUP_COLOR_CYCLE[i % len(GROUP_COLOR_CYCLE)]
                    for i, g in enumerate(group_order)}
    if len(group_order) > len(GROUP_COLOR_CYCLE):
        print(f"WARNING: {len(group_order)} group levels but only "
              f"{len(GROUP_COLOR_CYCLE)} default colors defined — colors will repeat.")

    # Layout: [row dendrogram | heatmap+strip | colorbar], with the middle
    # column split into a thin group-annotation strip above the heatmap
    # proper. The dendrogram column gets a matching blank top row so its
    # dendrogram aligns vertically with the heatmap (not the strip).
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
    ax_dend.set_ylim(ax_dend.get_ylim())  # lock before hiding ticks
    ax_dend.axis("off")

    # --- group annotation strip (top) ---
    # Drawn as individual axvspan rectangles (vector, one patch per group
    # block) rather than imshow — simpler than pcolormesh here since there
    # are only len(group_order) blocks, not one cell per sample.
    ax_strip.set_xlim(0, len(col_order))
    ax_strip.set_ylim(0, 1)
    ax_strip.set_xticks([]); ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)
    start = 0
    for g in group_order:
        n = int((groups.loc[col_order] == g).sum())
        ax_strip.axvspan(start, start + n, color=group_colors[g])
        ax_strip.text(start + n / 2, 0.5, f"{g} (n={n})", ha="center", va="center",
                      fontsize=6, fontweight="bold", color="white")
        start += n

    # --- heatmap ---
    # pcolormesh (not imshow) so every sample x pathway cell is drawn as an
    # individual vector polygon rather than a rasterized bitmap — stays crisp
    # at any zoom level in the saved PDF, and each cell is a distinct,
    # selectable object (e.g. editable in Illustrator). Thin white
    # edgecolors give each cell a visible border ("show each individual box
    # for each sample"). pcolormesh places row 0 at the BOTTOM by default
    # (opposite of imshow's default top-to-bottom convention) — invert_yaxis()
    # restores the same visual order as before (row_order top-to-bottom,
    # matching the y-tick labels).
    vlim = min(3.0, np.nanmax(np.abs(ordered.values))) if ordered.size else 1.0
    im = ax_heat.pcolormesh(ordered.values, cmap="RdBu_r", vmin=-vlim, vmax=vlim,
                            edgecolors="white", linewidth=0.4)
    ax_heat.invert_yaxis()
    # Cell (row i, col j) spans x in [j, j+1] under pcolormesh's default grid,
    # so a boundary BETWEEN column b-1 and column b sits at x=b exactly (no
    # -0.5 offset needed here, unlike imshow's pixel-center convention).
    for b in col_boundaries:
        ax_heat.axvline(b, color="black", linewidth=2.2)
    ax_heat.set_xticks([])
    ax_heat.set_yticks(np.arange(len(row_order)) + 0.5)
    ax_heat.set_yticklabels(
        [wrap_feature_label(f, width=label_wrap_width) for f in row_order], fontsize=6)
    ax_heat.set_xlabel(f"Samples (n={len(col_order)}, clustered within group)", fontsize=7)

    # --- colorbar (own gridspec cell — avoids manual cbar_pos overlap issues) ---
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Row z-score\n(log$_{10}$ relative abundance)", fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    ax_heat.set_title("Pathway enrichment across samples\n"
                      "(rows clustered; columns split by group, clustered within group)",
                      loc="left", fontweight="bold", fontsize=8)


# -----------------------------------------------------------------------------
# Main figure assembly
# -----------------------------------------------------------------------------
def make_figure(feat_df: pd.DataFrame, metadata: str, qval_threshold: float,
                top_n: int, label_top_n: int, output_prefix: str,
                label_wrap_width: int = 28,
                z_mat: pd.DataFrame = None, groups: pd.Series = None,
                group_order: list = None):
    set_publication_style()

    top = feat_df.nsmallest(top_n, "qval_isolated").iloc[::-1]  # smallest q at top of plot

    include_heatmap = z_mat is not None and groups is not None
    if include_heatmap:
        fig = plt.figure(figsize=(7.6, 9.0), constrained_layout=True)
        gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.5])
        top_row = gs[0].subgridspec(1, 2, wspace=0.35)
        ax_volcano = fig.add_subplot(top_row[0])
        ax_bar = fig.add_subplot(top_row[1])
    else:
        fig = plt.figure(figsize=(7.2, 3.4), constrained_layout=True)
        gs = fig.add_gridspec(1, 2)
        ax_volcano = fig.add_subplot(gs[0])
        ax_bar = fig.add_subplot(gs[1])

    plot_volcano(ax_volcano, feat_df, qval_threshold, label_top_n)
    plot_top_features(ax_bar, top, top_n, qval_threshold, label_wrap_width=label_wrap_width)

    if include_heatmap:
        # Panel C uses the SAME top pathways as Panel B, for direct
        # cross-reference between the bar chart and the heatmap rows.
        top_features_for_heatmap = [f for f in top["feature"] if f in z_mat.index]
        missing = set(top["feature"]) - set(top_features_for_heatmap)
        if missing:
            print(f"WARNING: {len(missing)} top-panel-B feature(s) not present in "
                 f"the abundance table and omitted from Panel C: {sorted(missing)}")
        plot_heatmap(fig, gs[1], z_mat.loc[top_features_for_heatmap], groups,
                    group_order, label_wrap_width=label_wrap_width)

    fig.suptitle(f"MaAsLin3 associations with {metadata}", fontsize=9, y=1.01)

    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Path to MaAsLin3 all_results.tsv")
    ap.add_argument("--metadata", required=True,
                    help="Metadata term to isolate, e.g. sarc_status_bin")
    ap.add_argument("--qval-threshold", type=float, default=0.1,
                    help="Significance line drawn on the volcano plot (default 0.1, matching MAASLIN_MAX_SIGNIFICANCE)")
    ap.add_argument("--top-n", type=int, default=15,
                    help="Number of features shown in the bar panel and heatmap (default 15)")
    ap.add_argument("--label-top-n", type=int, default=8,
                    help="Number of points labeled directly on the volcano plot (default 8)")
    ap.add_argument("--output", default="maaslin_volcano_top_features",
                    help="Output file prefix (default maaslin_volcano_top_features)")
    ap.add_argument("--omit-pathways", default="",
                    help="Comma-separated pathway IDs or full feature names to "
                         "exclude from the table and figure, e.g. "
                         "'PWY-101,PWY-6737: starch degradation V'")
    ap.add_argument("--omit-file", default=None,
                    help="Path to a text file with one pathway ID or full "
                         "feature name to omit per line. Combined with "
                         "--omit-pathways if both are given.")
    ap.add_argument("--label-wrap-width", type=int, default=28,
                    help="Character width at which pathway labels wrap onto "
                         "a new line, rather than being truncated (default 28)")
    ap.add_argument("--abundance-table", default=None,
                    help="Filtered feature abundance table (feature x sample; "
                         "see module docstring for the R export snippet). "
                         "Required, along with --sample-metadata, to draw Panel C.")
    ap.add_argument("--sample-metadata", default=None,
                    help="Sample -> group metadata table. Required, along with "
                         "--abundance-table, to draw Panel C.")
    ap.add_argument("--sample-id-col", default="File_ID",
                    help="Sample identifier column name in --sample-metadata (default File_ID)")
    ap.add_argument("--group-col", default=None,
                    help="Group column name in --sample-metadata used to split/color "
                         "heatmap columns (default: same as --metadata)")
    ap.add_argument("--group-order", default=None,
                    help="Comma-separated group level order for the heatmap columns "
                         "(default: sorted unique values observed)")
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    required_cols = {"feature", "metadata", "coef", "pval_individual", "pval_joint",
                     "error", "model", "N_not_zero"}
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(f"Input file is missing expected column(s): {missing}")

    omit_list = [s.strip() for s in args.omit_pathways.split(",") if s.strip()]
    if args.omit_file:
        with open(args.omit_file) as f:
            omit_list += [line.strip() for line in f if line.strip()]

    feat_df = prepare_feature_table(df, args.metadata)
    feat_df = omit_pathways(feat_df, omit_list)
    feat_df.to_csv(f"{args.output}_table.csv", index=False)

    z_mat = groups = group_order = None
    if args.abundance_table or args.sample_metadata:
        if not (args.abundance_table and args.sample_metadata):
            sys.exit("--abundance-table and --sample-metadata must both be given "
                     "to draw Panel C (or neither, to skip it).")
        group_col = args.group_col or args.metadata
        abundance_df = load_abundance_table(args.abundance_table)
        meta_df = load_sample_metadata(args.sample_metadata, args.sample_id_col, group_col)

        top_for_heatmap = feat_df.nsmallest(args.top_n, "qval_isolated")["feature"].tolist()
        z_mat, groups = build_enrichment_matrix(abundance_df, meta_df, top_for_heatmap, group_col)

        if args.group_order:
            group_order = [s.strip() for s in args.group_order.split(",")]
            unknown = set(group_order) - set(groups.unique())
            if unknown:
                sys.exit(f"--group-order contains level(s) not present in the data: {unknown}")
        else:
            group_order = sorted(groups.unique().tolist())

    make_figure(feat_df, args.metadata, args.qval_threshold,
               args.top_n, args.label_top_n, args.output,
               label_wrap_width=args.label_wrap_width,
               z_mat=z_mat, groups=groups, group_order=group_order)

    n_sig = (feat_df["qval_isolated"] < args.qval_threshold).sum()
    print(f"{len(feat_df)} features plotted. "
         f"{n_sig} significant at q < {args.qval_threshold} "
         f"(isolated to '{args.metadata}' only, not the full multi-covariate correction).")
    panel_c_note = " + heatmap (Panel C)" if z_mat is not None else " (Panel C skipped — no --abundance-table/--sample-metadata given)"
    print(f"Wrote {args.output}.pdf, {args.output}.png, {args.output}_table.csv{panel_c_note}")


if __name__ == "__main__":
    main()