#!/usr/bin/env python3
"""
plot_maaslin_results.py

Two-panel figure summarizing a MaAsLin3 all_results.tsv
for a single metadata term of interest (e.g. sarc_status_bin):

  Panel A: volcano plot (abundance-model coefficient vs. -log10(q-value)),
           colored by which sub-model (abundance vs. prevalence) drove the
           joint significance for that feature.
  Panel B: horizontal bar chart of the top-N features by q-value, signed by
           effect direction, colored the same way as Panel A for a
           consistent visual language across the figure.

Why q-values are recomputed here rather than trusting MaAsLin3's own
qval_joint column directly: qval_joint is FDR-corrected across every
metadata term tested (all fixed effects, not just the one of interest),
which makes it needlessly conservative for a single-variable figure. This
script isolates the metadata term with --metadata, then BH-corrects just
that term's pval_joint across features — matching the guidance in
MaAsLin3's own documentation ("it may be preferable to FDR correct just the
p-values from the variables of interest").

Usage:
    python plot_maaslin_results.py \
        --input all_results.tsv \
        --metadata sarc_status_bin \
        --qval-threshold 0.1 \
        --top-n 15 \
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
from matplotlib.lines import Line2D


# -----------------------------------------------------------------------------
# BH (Benjamini-Hochberg) correction, implemented directly rather than via
# statsmodels/scipy, so this script has no dependency beyond
# pandas/numpy/matplotlib. Matches R's p.adjust(method = "BH").
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


DRIVER_COLORS = {"abundance": "#2166AC", "prevalence": "#B2182B"}
DRIVER_LABELS = {"abundance": "Abundance-driven", "prevalence": "Prevalence-driven"}


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
    # Anchored in axes-fraction x (via get_yaxis_transform: x in [0,1] relative
    # to the axes, y in data coordinates) so the label always stays inside the
    # panel regardless of axis limits, instead of extending past the plot in
    # data coordinates and forcing constrained_layout to add huge inter-panel
    # spacing to avoid overlapping panel B.
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
def plot_top_features(ax, feat_df: pd.DataFrame, top_n: int, qval_threshold: float):
    top = feat_df.nsmallest(top_n, "qval_isolated").iloc[::-1]  # smallest q at top of plot
    labels = [clean_feature_label(f) for f in top["feature"]]
    colors = [DRIVER_COLORS[d] for d in top["driver"]]
    values = top["abund_coef"].values
    ypos = np.arange(len(top))

    ax.barh(ypos, values, color=colors, height=0.65, edgecolor="white", linewidth=0.4)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.axvline(0, color="grey", linewidth=0.6)

    for y, (_, row) in zip(ypos, top.iterrows()):
        sig_marker = "*" if row["qval_isolated"] < qval_threshold else ""
        x_text = row["abund_coef"] + (0.02 if row["abund_coef"] >= 0 else -0.02)
        ha = "left" if row["abund_coef"] >= 0 else "right"
        ax.text(x_text, y, f"q={row['qval_isolated']:.2g}{sig_marker}",
                va="center", ha=ha, fontsize=6, color="black")

    ax.set_xlabel("Abundance-model coefficient (log$_2$ fold change)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"Top {top_n} pathways by q-value", loc="left", fontweight="bold", fontsize=8)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def make_figure(feat_df: pd.DataFrame, metadata: str, qval_threshold: float,
                top_n: int, label_top_n: int, output_prefix: str):
    set_publication_style()

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)

    plot_volcano(axes[0], feat_df, qval_threshold, label_top_n)
    plot_top_features(axes[1], feat_df, top_n, qval_threshold)

    fig.suptitle(f"MaAsLin3 associations with {metadata}", fontsize=9, y=1.06)

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
                    help="Number of features shown in the bar panel (default 15)")
    ap.add_argument("--label-top-n", type=int, default=8,
                    help="Number of points labeled directly on the volcano plot (default 8)")
    ap.add_argument("--output", default="maaslin_volcano_top_features",
                    help="Output file prefix (default maaslin_volcano_top_features)")
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    required_cols = {"feature", "metadata", "coef", "pval_individual", "pval_joint",
                     "error", "model", "N_not_zero"}
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(f"Input file is missing expected column(s): {missing}")

    feat_df = prepare_feature_table(df, args.metadata)
    feat_df.to_csv(f"{args.output}_table.csv", index=False)

    make_figure(feat_df, args.metadata, args.qval_threshold,
               args.top_n, args.label_top_n, args.output)

    n_sig = (feat_df["qval_isolated"] < args.qval_threshold).sum()
    print(f"{len(feat_df)} features plotted. "
         f"{n_sig} significant at q < {args.qval_threshold} "
         f"(isolated to '{args.metadata}' only, not the full multi-covariate correction).")
    print(f"Wrote {args.output}.pdf, {args.output}.png, {args.output}_table.csv")


if __name__ == "__main__":
    main()