"""
plot_significance_grid.py

Publication-quality figure: features (rows) x models (columns), with a
colored box wherever that model calls the feature significant. For DESeq2
and lasso logistic regression, box color also encodes the direction of the
effect (increased vs decreased in cases). A count of significant features
is printed beneath each model's column.

Significance rules (as specified):
  DESeq2   : deseq2_padj < 0.05
             direction: sign(deseq2_log2FC)
  Lasso    : exp(linreg_lower_perc12_5), exp(linreg_upper_perc87_5) -> odds
             range does NOT include 1
             direction: sign(linreg_mean_coef)
  FCNN     : [fcnn_lower_perc12_5, fcnn_upper_perc87_5] does NOT include 0
  RF       : rf_selection_frequency > 0

By default, only features significant in at least 2 of the 4 models are
plotted (--min-sig to change). The per-column counts shown beneath the grid
are the TOTAL number of significant features for that model (i.e. before
the min-sig filter is applied), so they reflect each model's overall yield.

Input: one merged CSV produced by feature_importance.merge_importance_tables(),
containing at least the columns:
  feature, deseq2_padj, deseq2_log2FC,
  linreg_mean_coef, linreg_lower_perc12_5, linreg_upper_perc87_5,
  fcnn_lower_perc12_5, fcnn_upper_perc87_5,
  rf_selection_frequency

Usage:
  python plot_significance_grid.py <merged_csv> [--out OUT_PREFIX]
                                    [--min-sig N] [--sort-by count]
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

MODEL_COLS = {
    "DESeq2": "sig_deseq2",
    "Lasso logistic regression": "sig_lasso",
    "FCNN": "sig_fcnn",
    "Random forest": "sig_rf",
}

# Models whose box color encodes effect direction vs. models with a single
# "significant" color (no direction requested for these).
DIRECTIONAL_MODELS = {"DESeq2", "Lasso logistic regression"}

POS_COLOR = "#B23A48"   # increased in cases
NEG_COLOR = "#2A6F97"   # decreased in cases
FCNN_COLOR = "#55A868"
RF_COLOR = "#8E6C9E"
NONSIG_FACE = "#F2F2F2"
NONSIG_EDGE = "#D9D9D9"

LEGEND_ITEMS = [
    Patch(facecolor=POS_COLOR, edgecolor=POS_COLOR, label="Increased in cases (DESeq2 / Lasso)"),
    Patch(facecolor=NEG_COLOR, edgecolor=NEG_COLOR, label="Decreased in cases (DESeq2 / Lasso)"),
    Patch(facecolor=FCNN_COLOR, edgecolor=FCNN_COLOR, label="Significant (FCNN)"),
    Patch(facecolor=RF_COLOR, edgecolor=RF_COLOR, label="Significant (Random forest)"),
    Patch(facecolor=NONSIG_FACE, edgecolor=NONSIG_EDGE, label="Not significant"),
]


def compute_significance(df):
    df = df.copy()

    df["sig_deseq2"] = df["deseq2_padj"] < 0.05
    df["dir_deseq2"] = np.sign(df["deseq2_log2FC"])

    lo_or = np.exp(df["linreg_lower_perc12_5"])
    hi_or = np.exp(df["linreg_upper_perc87_5"])
    df["sig_lasso"] = ~((lo_or <= 1) & (hi_or >= 1))
    df["dir_lasso"] = np.sign(df["linreg_mean_coef"])

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


def cell_color(row, model):
    sig_col = MODEL_COLS[model]
    if not row[sig_col]:
        return NONSIG_FACE, NONSIG_EDGE

    if model == "DESeq2":
        color = POS_COLOR if row["dir_deseq2"] > 0 else NEG_COLOR
    elif model == "Lasso logistic regression":
        color = POS_COLOR if row["dir_lasso"] > 0 else NEG_COLOR
    elif model == "FCNN":
        color = FCNN_COLOR
    else:  # Random forest
        color = RF_COLOR
    return color, color


def plot_grid(df, out_prefix="significance_grid", min_sig=2, sort_by="count"):
    df = compute_significance(df)

    # Total significant-feature counts per model, computed BEFORE the
    # min-sig filter, so they reflect each model's full yield.
    total_counts = {model: int(df[col].sum()) for model, col in MODEL_COLS.items()}

    df = df[df["n_sig"] >= min_sig].copy()

    if sort_by == "count":
        df = df.sort_values(["n_sig", "feature"], ascending=[False, True])
    else:
        df = df.sort_values("feature")

    if len(df) == 0:
        raise ValueError(
            f"No features significant in at least {min_sig} models.")

    df["feature_label"] = df["feature"].apply(clean_feature_name)

    models = list(MODEL_COLS.keys())
    n_feat = len(df)
    n_model = len(models)

    fig_h = max(3, 0.28 * n_feat + 2.2)
    fig_w = 2.2 * n_model + 3.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for row_i, (_, row) in enumerate(df.iterrows()):
        y = n_feat - row_i - 1
        for col_i, model in enumerate(models):
            face, edge = cell_color(row, model)
            ax.add_patch(Rectangle((col_i, y), 0.9, 0.9,
                                    facecolor=face, edgecolor=edge,
                                    linewidth=0.8))

    ax.set_xlim(0, n_model)
    ax.set_ylim(-1, n_feat)
    ax.set_xticks([i + 0.45 for i in range(n_model)])
    ax.set_xticklabels(models, fontsize=10, fontweight="bold",
                        rotation=30, ha="left", rotation_mode="anchor")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_yticks([n_feat - i - 0.55 for i in range(n_feat)])
    ax.set_yticklabels(df["feature_label"], fontsize=8, fontstyle="italic")

    # Count of total significant features beneath each model's column.
    for col_i, model in enumerate(models):
        ax.text(col_i + 0.45, -0.5, f"n = {total_counts[model]}",
                ha="center", va="center", fontsize=9)

    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    plt.title(f"Features significant in \u2265{min_sig} models",
              fontsize=13, fontweight="bold", pad=14)
    ax.legend(handles=LEGEND_ITEMS, loc="upper center",
              bbox_to_anchor=(0.5, -0.06 - 0.015 * n_feat if n_feat > 20 else -0.08),
              ncol=2, frameon=False, fontsize=8)
    plt.tight_layout()

    fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches="tight")
    print(f"Saved {out_prefix}.pdf and {out_prefix}.png "
          f"({n_feat} features, min_sig={min_sig})")
    print("Total significant features per model:", total_counts)
    return fig


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("merged_csv")
    ap.add_argument("--out", default="significance_grid")
    ap.add_argument("--min-sig", type=int, default=2,
                     help="Minimum number of models a feature must be "
                          "significant in to be plotted (default: 2).")
    ap.add_argument("--sort-by", choices=["count", "alpha"], default="count")
    args = ap.parse_args()

    df = pd.read_csv(args.merged_csv)
    plot_grid(df, out_prefix=args.out, min_sig=args.min_sig,
              sort_by=args.sort_by)