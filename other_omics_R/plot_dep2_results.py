#!/usr/bin/env python3
"""
plot_dep2_results.py

Three-panel figure summarizing a DEP2 differential-abundance results table
(the "<experiment>_dep2_results.csv" produced by 02_proteomics_dep2.R) for a
single contrast (e.g. Sarc_vs_NoSarc), across any of the DEP2-based assays:
GC-MS, LC-MS (pos/neg), or the massspec (DIA proteomics) loader.

  Panel A: volcano plot (diff vs. -log10(p-value)), colored by whether the
           feature required imputation (--> "imputed" column from DEP2),
           since that's the closest DEP2-native reliability signal to flag
           on a per-feature basis.
  Panel B: horizontal bar chart of the top-N features by q-value, signed by
           effect direction, with CI.L/CI.R error bars (DEP2 reports these
           natively; unlike the MaAsLin3 version there's no need to omit
           them), colored the same way as Panel A, with p- and q-values
           reported to the right of each bar.
  Panel C (optional): clustered heatmap of per-sample values (row z-scored)
           for the same top-N features shown in Panel B. Features (rows) are
           hierarchically clustered; samples (columns) are first split into
           groups by --group-col (e.g. sarc_status_bin), then hierarchically
           clustered *within* each group, so the group split is always
           visually respected and clustering only reorders samples inside a
           group.

Panel C requires two additional inputs beyond the results CSV, both already
produced by 02_proteomics_dep2.R without any extra export step:
  --abundance-table  <experiment>_dep2_vsn_imputed_matrix.csv   (sample x feature)
  --sample-metadata  <experiment>_dep2_vsn_imputed_metadata.csv (label + covariates)
If neither is supplied, the script falls back to the two-panel (A+B) figure.

Why q-values are recomputed here rather than trusting DEP2's own p.adj
column directly: DEP2's p.adj is BH-corrected over the full feature set that
went into test_diff(). If --omit-features/--omit-file drops any rows before
plotting, the correct multiple-testing correction is over the REMAINING
feature set, not the original one — so this script always recomputes BH
from the raw p-value column after any omissions, and keeps DEP2's own
p.adj only as a reference column in the output table (for auditing/sanity
checking, printed as a discrepancy warning if they diverge).

Usage:
    python plot_dep2_results.py \
        --input GC_MS_dep2_results.csv \
        --contrast Sarc_vs_NoSarc \
        --experiment-name GC-MS \
        --qval-threshold 0.1 \
        --top-n 15 \
        --abundance-table GC_MS_dep2_vsn_imputed_matrix.csv \
        --sample-metadata GC_MS_dep2_vsn_imputed_metadata.csv \
        --sample-id-col label \
        --group-col sarc_status_bin \
        --output GC_MS_dep2_volcano_top_features

Produces <output>.pdf and <output>.png (300 DPI), plus <output>_table.csv
with the recomputed q-values, so the numbers behind the figure are
auditable rather than baked silently into a plot.
"""

import argparse
import re
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
# pandas/numpy/matplotlib/scipy. Matches R's p.adjust(method = "BH") — and
# DEP2's own fdr.type = "BH" in test_diff(), so recomputing it here (see
# module docstring) reproduces DEP2's p.adj exactly when no features have
# been omitted.
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
# Load + reshape a DEP2 results CSV down to one row per feature for the
# contrast of interest. Unlike MaAsLin3's all_results.tsv, DEP2's output is
# already one row per feature with a single test per contrast (no
# abundance/prevalence submodel split), so this is a column-selection +
# rename rather than a long-to-wide reshape.
# -----------------------------------------------------------------------------
def prepare_feature_table(df: pd.DataFrame, contrast: str) -> pd.DataFrame:
    diff_col = f"{contrast}_diff"
    pval_col = f"{contrast}_p.val"
    padj_col = f"{contrast}_p.adj"
    ci_l_col = f"{contrast}_CI.L"
    ci_r_col = f"{contrast}_CI.R"

    required = {diff_col, pval_col, ci_l_col, ci_r_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Input is missing expected column(s) for contrast '{contrast}': "
            f"{sorted(missing)}. Available columns: {list(df.columns)}"
        )

    feature_col = "name" if "name" in df.columns else (
        "ID" if "ID" in df.columns else df.columns[0])

    out = pd.DataFrame({
        "feature": df[feature_col].astype(str),
        "diff": df[diff_col].astype(float),
        "pval": df[pval_col].astype(float),
        "ci_l": df[ci_l_col].astype(float),
        "ci_r": df[ci_r_col].astype(float),
    })
    out["imputed"] = df["imputed"].astype(bool) if "imputed" in df.columns else False
    out["num_nas"] = df["num_NAs"] if "num_NAs" in df.columns else np.nan
    if padj_col in df.columns:
        # DEP2's own p.adj — kept only as a reference/sanity-check column,
        # not used for plotting. See module docstring.
        out["padj_dep2"] = df[padj_col].astype(float)

    n_before = len(out)
    out = out.dropna(subset=["pval"]).reset_index(drop=True)
    if len(out) < n_before:
        print(f"NOTE: dropped {n_before - len(out)} feature(s) with missing "
              f"{pval_col} before plotting.")

    return out.sort_values("pval").reset_index(drop=True)


def finalize_qvalues(feat_df: pd.DataFrame) -> pd.DataFrame:
    """Recompute BH q-values (see module docstring) and, if DEP2's own
    p.adj is present, warn if it disagrees non-trivially with the recompute
    — expected only when features were omitted (changes n) or never
    (sanity check that the manual BH implementation matches DEP2's)."""
    feat_df = feat_df.copy()
    feat_df["qval"] = bh_qvalues(feat_df["pval"].values)
    if "padj_dep2" in feat_df.columns:
        diff = (feat_df["qval"] - feat_df["padj_dep2"]).abs()
        n_diverge = int((diff > 1e-6).sum())
        if n_diverge > 0:
            print(f"NOTE: recomputed q-values differ from DEP2's own p.adj "
                  f"for {n_diverge}/{len(feat_df)} feature(s) — expected if "
                  f"--omit-features/--omit-file dropped any rows (changes "
                  f"the multiple-testing set), unexpected otherwise.")
    return feat_df.sort_values("pval").reset_index(drop=True)


# -----------------------------------------------------------------------------
# omit_features — drop rows whose feature exactly matches an entry in
# omit_list. Exact match only (not the MaAsLin3 version's colon-prefix
# matching, since DEP2 feature names don't follow a "ID: name" convention
# consistently across assays — safer to require an exact copy-paste).
# -----------------------------------------------------------------------------
def omit_features(df: pd.DataFrame, omit_list: list) -> pd.DataFrame:
    if not omit_list:
        return df
    omit_set = set(omit_list)
    keep_mask = ~df["feature"].isin(omit_set)
    n_dropped = (~keep_mask).sum()
    if n_dropped > 0:
        print(f"Omitting {n_dropped} feature(s) per --omit-features/--omit-file: "
              f"{df.loc[~keep_mask, 'feature'].tolist()}")
    unmatched = omit_set - set(df["feature"])
    if unmatched:
        print(f"WARNING: {len(unmatched)} omit entry(ies) did not match any "
              f"feature in the input and had no effect: {sorted(unmatched)}")
    return df[keep_mask].reset_index(drop=True)


# -----------------------------------------------------------------------------
# Cosmetic helpers
# -----------------------------------------------------------------------------
# Matches the GC-MS/LC-MS wide-loader naming convention, e.g.
# "[10394] 3-(4-hydroxyphenyl)propionic acid [15.979]" -> the compound name
# in the middle (bracketed feature ID .. bracketed retention time).
_BRACKETED_ID_RT_RE = re.compile(r"^\[[^\]]+\]\s*(.*?)\s*\[[^\]]+\]$")


def clean_feature_label(feature: str, max_len: int = 45) -> str:
    """Produce a readable display label from a DEP2 feature identifier.
    Naming conventions vary by assay/loader:
      - GC-MS / LC-MS (wide loader): "[ID] compound name [RT]"
        -> strip the leading bracketed ID and trailing bracketed RT.
      - massspec loader: "ProteinGroup;GeneName"
        -> keep the gene name (text after the last ';').
      - LC-MS clustered features (03_lc_feature_exploration.R output):
        "Cluster12" -> used as-is.
    Anything that doesn't match one of these falls back to the raw feature
    string unchanged, so an unfamiliar naming convention degrades
    gracefully rather than mangling the label."""
    label = feature.strip()
    m = _BRACKETED_ID_RT_RE.match(label)
    if m and m.group(1):
        label = m.group(1).strip()
    elif ";" in label:
        label = label.rsplit(";", 1)[-1].strip()
    if len(label) > max_len:
        label = label[: max_len - 1] + "\u2026"
    return label


def wrap_feature_label(feature: str, width: int = 28) -> str:
    """Like clean_feature_label, but wraps onto multiple lines (via
    textwrap.fill) instead of truncating with an ellipsis — for contexts
    (e.g. the bar chart's y-axis) where the full name is preferred and
    vertical space for extra lines is available."""
    label = feature.strip()
    m = _BRACKETED_ID_RT_RE.match(label)
    if m and m.group(1):
        label = m.group(1).strip()
    elif ";" in label:
        label = label.rsplit(";", 1)[-1].strip()
    return textwrap.fill(label, width=width)


IMPUTED_COLORS = {False: "#2166AC", True: "#B2182B"}
IMPUTED_LABELS = {False: "Fully observed", True: "Contains imputed value(s)"}

# Default group color cycle for Panel C's group-annotation strip. Extended
# as needed if --group-col has more than the number of colors below (falls
# back to a warning and repeats colors).
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


def diff_axis_label(contrast: str) -> str:
    """'Sarc_vs_NoSarc' -> 'Sarc \u2212 NoSarc (VSN-normalized log2 intensity
    diff.)' — falls back to a generic label if the contrast name doesn't
    follow DEP2's '<A>_vs_<B>' test_diff() naming convention."""
    m = re.match(r"^(.+)_vs_(.+)$", contrast)
    if m:
        return f"{m.group(1)} \u2212 {m.group(2)} (VSN-normalized log$_2$ intensity diff.)"
    return f"{contrast} diff. (VSN-normalized log$_2$ intensity)"


# -----------------------------------------------------------------------------
# Panel A: volcano plot
# -----------------------------------------------------------------------------
def compute_bh_pvalue_threshold(feat_df: pd.DataFrame, qval_threshold: float):
    """Convert a q-value significance threshold into a single p-value line for
    a raw-p volcano plot. Two cases, since BH significance is rank-dependent
    rather than a fixed p-value cutoff:

      - If >=1 feature already clears qval_threshold: the true empirical BH
        cutoff is the LARGEST p-value among those significant features —
        return that, with is_empirical=True.
      - If none do: there's no empirical cutoff to report. Fall back to the
        rank-1 BH critical value (alpha / n) — the smallest, strictest
        boundary, representing "how small the single best p-value would have
        needed to be." Returned with is_empirical=False so the plot can label
        it as a conservative reference rather than an observed cutoff.
    """
    n = len(feat_df)
    sig = feat_df[feat_df["qval"] < qval_threshold]
    if len(sig) > 0:
        return sig["pval"].max(), True
    return qval_threshold / n, False


def plot_volcano(ax, feat_df: pd.DataFrame, qval_threshold: float, label_top_n: int,
                 contrast: str):
    x = feat_df["diff"].values
    y = -np.log10(feat_df["pval"].values.clip(min=1e-300))

    for imputed, sub in feat_df.groupby("imputed"):
        idx = feat_df.index.isin(sub.index)
        ax.scatter(
            x[idx], y[idx],
            s=22, linewidths=0.4, edgecolors="white",
            color=IMPUTED_COLORS[imputed], alpha=0.85,
            label=IMPUTED_LABELS[imputed], zorder=3,
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

    top = feat_df.nsmallest(label_top_n, "qval")
    for _, row in top.iterrows():
        ax.annotate(
            clean_feature_label(row["feature"], max_len=28),
            xy=(row["diff"], -np.log10(max(row["pval"], 1e-300))),
            xytext=(4, 3), textcoords="offset points",
            fontsize=6, color="black",
        )

    ax.set_xlabel(diff_axis_label(contrast))
    ax.set_ylabel(r"$-\log_{10}(p\mathrm{-value})$")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper left", handletextpad=0.4, borderaxespad=0.2)


# -----------------------------------------------------------------------------
# Panel B: top-N features, horizontal bar chart
# -----------------------------------------------------------------------------
def plot_top_features(ax, top: pd.DataFrame, top_n: int, qval_threshold: float,
                      contrast: str, label_wrap_width: int = 28):
    labels = [wrap_feature_label(f, width=label_wrap_width) for f in top["feature"]]
    colors = [IMPUTED_COLORS[imp] for imp in top["imputed"]]
    values = top["diff"].values

    # Asymmetric CI error bars, unique to DEP2's per-feature CI.L/CI.R
    # (MaAsLin3's all_results.tsv doesn't carry these). Clip at 0 to guard
    # against tiny floating-point rounding putting ci_l fractionally above
    # diff or vice versa.
    err_lower = np.clip(values - top["ci_l"].values, 0, None)
    err_upper = np.clip(top["ci_r"].values - values, 0, None)

    max_lines = max(label.count("\n") + 1 for label in labels)
    row_spacing = 1.0 + 0.45 * (max_lines - 1)
    ypos = np.arange(len(top)) * row_spacing
    bar_height = 0.65 * min(row_spacing, 1.6)

    ax.barh(ypos, values, xerr=[err_lower, err_upper], color=colors, height=bar_height,
           edgecolor="white", linewidth=0.4,
           error_kw=dict(elinewidth=0.7, capsize=2, capthick=0.7, ecolor="#444444"))
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_ylim(-row_spacing * 0.75, ypos[-1] + row_spacing * 0.75)
    ax.axvline(0, color="grey", linewidth=0.6)

    for y, (_, row) in zip(ypos, top.iterrows()):
        sig_marker = "*" if row["qval"] < qval_threshold else ""
        ax.text(1.03, y, f"p={row['pval']:.2g}, q={row['qval']:.2g}{sig_marker}",
                transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=6, color="black")

    ax.set_xlabel(diff_axis_label(contrast))
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"Top {top_n} features by q-value", loc="left", fontweight="bold", fontsize=8)


# -----------------------------------------------------------------------------
# Panel C: clustered heatmap of per-sample values for the top features
# -----------------------------------------------------------------------------
def load_abundance_table(path: str) -> pd.DataFrame:
    """DEP2's <experiment>_dep2_vsn_imputed_matrix.csv is sample x feature
    (a 'Sample' column, then one column per feature) — the opposite
    orientation from the internal feature x sample convention used below,
    so this loader transposes after reading."""
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    sample_col = "Sample" if "Sample" in df.columns else df.columns[0]
    df = df.set_index(sample_col)
    df.index.name = "Sample"
    df.index = df.index.astype(str)
    return df.T  # -> feature x sample


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
    """Restrict to the top features x samples with known group membership,
    then row-z-score. Unlike the MaAsLin3 version, there's no
    log10(x + pseudocount) step here: DEP2's imputed matrix is already
    VSN-normalized (a log2-like scale) and fully imputed (no zeros/missing
    values by construction, since se_imp has already been through
    DEP2::impute()), so a second log transform on top would distort the
    scale rather than stabilize it."""
    abundance_df = abundance_df.copy()
    abundance_df.columns = abundance_df.columns.astype(str)

    missing_features = [f for f in top_features if f not in abundance_df.index]
    if missing_features:
        sys.exit(
            "The following top-panel-B feature(s) were not found in "
            "--abundance-table (row identifiers must match the results "
            f"CSV's 'name'/'ID' column exactly): {missing_features}"
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

    row_mean = mat.mean(axis=1)
    row_std = mat.std(axis=1, ddof=0)
    # Guard against a (rare, post-filtering) zero-variance row: z-score would
    # be 0/0. Leave it as a flat 0 (no enrichment signal) rather than NaN,
    # so it doesn't break clustering distances.
    row_std_safe = row_std.replace(0, np.nan)
    z_mat = mat.sub(row_mean, axis=0).div(row_std_safe, axis=0).fillna(0.0)

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
    # pcolormesh (not imshow) so every sample x feature cell is drawn as an
    # individual vector polygon — stays crisp at any zoom level in the saved
    # PDF, and each cell is a distinct, selectable object. pcolormesh places
    # row 0 at the BOTTOM by default (opposite of imshow) — invert_yaxis()
    # restores top-to-bottom order matching the y-tick labels.
    vlim = min(3.0, np.nanmax(np.abs(ordered.values))) if ordered.size else 1.0
    im = ax_heat.pcolormesh(ordered.values, cmap="RdBu_r", vmin=-vlim, vmax=vlim,
                            edgecolors="white", linewidth=0.4)
    ax_heat.invert_yaxis()
    for b in col_boundaries:
        ax_heat.axvline(b, color="black", linewidth=2.2)
    ax_heat.set_xticks([])
    ax_heat.set_yticks(np.arange(len(row_order)) + 0.5)
    ax_heat.set_yticklabels(
        [wrap_feature_label(f, width=label_wrap_width) for f in row_order], fontsize=6)
    ax_heat.set_xlabel(f"Samples (n={len(col_order)}, clustered within group)", fontsize=7)

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Row z-score\n(VSN-normalized log$_2$ intensity)", fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    ax_heat.set_title("Feature levels across samples\n"
                      "(rows clustered; columns split by group, clustered within group)",
                      loc="left", fontweight="bold", fontsize=8)


# -----------------------------------------------------------------------------
# Main figure assembly
# -----------------------------------------------------------------------------
def make_figure(feat_df: pd.DataFrame, contrast: str, experiment_name: str,
                qval_threshold: float, top_n: int, label_top_n: int, output_prefix: str,
                label_wrap_width: int = 28,
                z_mat: pd.DataFrame = None, groups: pd.Series = None,
                group_order: list = None):
    set_publication_style()

    top = feat_df.nsmallest(top_n, "qval").iloc[::-1]  # smallest q at top of plot

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

    plot_volcano(ax_volcano, feat_df, qval_threshold, label_top_n, contrast)
    plot_top_features(ax_bar, top, top_n, qval_threshold, contrast,
                      label_wrap_width=label_wrap_width)

    if include_heatmap:
        # Panel C uses the SAME top features as Panel B, for direct
        # cross-reference between the bar chart and the heatmap rows.
        top_features_for_heatmap = [f for f in top["feature"] if f in z_mat.index]
        missing = set(top["feature"]) - set(top_features_for_heatmap)
        if missing:
            print(f"WARNING: {len(missing)} top-panel-B feature(s) not present in "
                 f"the abundance table and omitted from Panel C: {sorted(missing)}")
        plot_heatmap(fig, gs[1], z_mat.loc[top_features_for_heatmap], groups,
                    group_order, label_wrap_width=label_wrap_width)

    title = f"DEP2 differential abundance: {contrast}"
    if experiment_name:
        title = f"{experiment_name} \u2014 {title}"
    fig.suptitle(title, fontsize=9, y=1.01)

    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True,
                    help="Path to <experiment>_dep2_results.csv")
    ap.add_argument("--contrast", default="Sarc_vs_NoSarc",
                    help="test_diff() contrast to plot, matching the "
                         "<contrast>_diff / _p.val / _CI.L / _CI.R column "
                         "prefix (default Sarc_vs_NoSarc)")
    ap.add_argument("--experiment-name", default=None,
                    help="Optional assay label for the figure title, e.g. 'GC-MS'")
    ap.add_argument("--qval-threshold", type=float, default=0.1,
                    help="Significance line drawn on the volcano plot (default 0.1)")
    ap.add_argument("--top-n", type=int, default=15,
                    help="Number of features shown in the bar panel and heatmap (default 15)")
    ap.add_argument("--label-top-n", type=int, default=8,
                    help="Number of points labeled directly on the volcano plot (default 8)")
    ap.add_argument("--output", default="dep2_volcano_top_features",
                    help="Output file prefix (default dep2_volcano_top_features)")
    ap.add_argument("--omit-features", default="",
                    help="Comma-separated exact feature names to exclude from "
                         "the table and figure")
    ap.add_argument("--omit-file", default=None,
                    help="Path to a text file with one exact feature name to "
                         "omit per line. Combined with --omit-features if both are given.")
    ap.add_argument("--label-wrap-width", type=int, default=28,
                    help="Character width at which feature labels wrap onto "
                         "a new line, rather than being truncated (default 28)")
    ap.add_argument("--abundance-table", default=None,
                    help="<experiment>_dep2_vsn_imputed_matrix.csv (sample x "
                         "feature). Required, along with --sample-metadata, "
                         "to draw Panel C.")
    ap.add_argument("--sample-metadata", default=None,
                    help="<experiment>_dep2_vsn_imputed_metadata.csv. "
                         "Required, along with --abundance-table, to draw Panel C.")
    ap.add_argument("--sample-id-col", default="label",
                    help="Sample identifier column name in --sample-metadata (default label)")
    ap.add_argument("--group-col", default="sarc_status_bin",
                    help="Group column name in --sample-metadata used to split/color "
                         "heatmap columns (default sarc_status_bin)")
    ap.add_argument("--group-order", default=None,
                    help="Comma-separated group level order for the heatmap columns "
                         "(default: sorted unique values observed)")
    args = ap.parse_args()

    df = pd.read_csv(args.input)

    omit_list = [s.strip() for s in args.omit_features.split(",") if s.strip()]
    if args.omit_file:
        with open(args.omit_file) as f:
            omit_list += [line.strip() for line in f if line.strip()]

    feat_df = prepare_feature_table(df, args.contrast)
    feat_df = omit_features(feat_df, omit_list)
    feat_df = finalize_qvalues(feat_df)
    feat_df.to_csv(f"{args.output}_table.csv", index=False)

    z_mat = groups = group_order = None
    if args.abundance_table or args.sample_metadata:
        if not (args.abundance_table and args.sample_metadata):
            sys.exit("--abundance-table and --sample-metadata must both be given "
                     "to draw Panel C (or neither, to skip it).")
        abundance_df = load_abundance_table(args.abundance_table)
        meta_df = load_sample_metadata(args.sample_metadata, args.sample_id_col, args.group_col)

        top_for_heatmap = feat_df.nsmallest(args.top_n, "qval")["feature"].tolist()
        z_mat, groups = build_enrichment_matrix(abundance_df, meta_df, top_for_heatmap,
                                                args.group_col)

        if args.group_order:
            group_order = [s.strip() for s in args.group_order.split(",")]
            unknown = set(group_order) - set(groups.unique())
            if unknown:
                sys.exit(f"--group-order contains level(s) not present in the data: {unknown}")
        else:
            group_order = sorted(groups.unique().tolist())

    make_figure(feat_df, args.contrast, args.experiment_name, args.qval_threshold,
               args.top_n, args.label_top_n, args.output,
               label_wrap_width=args.label_wrap_width,
               z_mat=z_mat, groups=groups, group_order=group_order)

    n_sig = (feat_df["qval"] < args.qval_threshold).sum()
    print(f"{len(feat_df)} features plotted. "
         f"{n_sig} significant at q < {args.qval_threshold}.")
    panel_c_note = " + heatmap (Panel C)" if z_mat is not None else " (Panel C skipped — no --abundance-table/--sample-metadata given)"
    print(f"Wrote {args.output}.pdf, {args.output}.png, {args.output}_table.csv{panel_c_note}")


if __name__ == "__main__":
    main()