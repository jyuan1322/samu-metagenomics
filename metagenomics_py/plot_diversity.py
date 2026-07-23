"""
plot_diversity.py

Two-panel publication figure:
  (A) Alpha diversity (Shannon index by default) per sample, as a boxplot
      split by Sarc / NonSarc, with individual points overlaid and a
      Mann-Whitney U p-value annotated.
  (B) Beta diversity via PCoA (classical/metric MDS, Torgerson's method) on
      Bray-Curtis dissimilarity between samples, colored by group, with 95%
      confidence ellipses per group and a PERMANOVA R^2 / p-value annotated
      (permutation test on the pseudo-F statistic, Anderson 2001).

Uses the SAME sample filtering as the rest of the pipeline (Full.SaMu==1,
age_def>=50, complete covariates) via loaders._clinical_covariates, but
operates on raw (untransformed) MetaPhlAn relative abundances -- diversity
metrics should be computed on proportions, not CLR- or log-transformed data.

Usage:
  python plot_diversity.py \
      --input-dir /data/local/jy1008/SaMu/results/latest/metagenomics_R \
      --meta-filtered-csv meta_filtered_06242026.csv \
      --meta-df-csv meta_df_FullSaMu_06242026.csv \
      --alpha-metric shannon \
      --n-permutations 999 \
      --out diversity_panels
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu

from loaders import _clinical_covariates

GROUP_COLORS = {"NonSarc": "#4C72B0", "Sarc": "#C44E52"}
GROUP_ORDER = ["NonSarc", "Sarc"]


# ---------------------------------------------------------------------------
# Alpha diversity metrics (computed per sample from a relative-abundance row)
# ---------------------------------------------------------------------------
def shannon(p):
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def simpson(p):
    return float(1 - (p ** 2).sum())


def richness(p):
    return float((p > 0).sum())


ALPHA_METRICS = {"shannon": shannon, "simpson": simpson, "richness": richness}


# ---------------------------------------------------------------------------
# Classical (metric) PCoA -- Torgerson's method
# ---------------------------------------------------------------------------
def pcoa(dist_matrix, n_axes=2):
    n = dist_matrix.shape[0]
    D2 = dist_matrix ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    total_pos = eigvals[eigvals > 0].sum()
    var_explained = np.divide(eigvals, total_pos, out=np.zeros_like(eigvals),
                              where=total_pos > 0)

    coords = eigvecs[:, :n_axes] * np.sqrt(np.clip(eigvals[:n_axes], 0, None))
    return coords, var_explained[:n_axes]


# ---------------------------------------------------------------------------
# PERMANOVA (Anderson 2001 pseudo-F, permutation p-value)
# ---------------------------------------------------------------------------
def permanova(dist_matrix, groups, n_perm=999, seed=0):
    groups = np.asarray(groups)
    labels = np.unique(groups)
    n = dist_matrix.shape[0]
    D2 = dist_matrix ** 2

    def pseudo_f(assign):
        ss_total = D2[np.triu_indices(n, k=1)].sum() / n
        ss_within = 0.0
        for lab in labels:
            idx = np.where(assign == lab)[0]
            if len(idx) < 2:
                continue
            sub = D2[np.ix_(idx, idx)]
            ss_within += sub[np.triu_indices(len(idx), k=1)].sum() / len(idx)
        ss_among = ss_total - ss_within
        a = len(labels)
        df_among, df_within = a - 1, n - a
        f_stat = (ss_among / df_among) / (ss_within / df_within)
        r2 = ss_among / ss_total
        return f_stat, r2

    f_obs, r2_obs = pseudo_f(groups)

    rng = np.random.default_rng(seed)
    count_ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(groups)
        f_perm, _ = pseudo_f(perm)
        if f_perm >= f_obs:
            count_ge += 1
    p_value = (count_ge + 1) / (n_perm + 1)
    return f_obs, r2_obs, p_value


def confidence_ellipse(ax, x, y, color, n_std=1.96, **kwargs):
    """Approximate 95% confidence ellipse (n_std=1.96) for a 2D point cloud,
    via eigendecomposition of the covariance matrix."""
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    width, height = 2 * n_std * np.sqrt(np.clip(eigvals, 0, None))
    ell = Ellipse((np.mean(x), np.mean(y)), width, height, angle=angle,
                 facecolor=color, alpha=0.12, edgecolor=color, linewidth=1.2,
                 **kwargs)
    ax.add_patch(ell)


def load_abundance(input_dir, meta_filtered_csv, meta_df_csv):
    input_dir = Path(input_dir)
    meta_filtered = pd.read_csv(input_dir / meta_filtered_csv)
    meta_df = pd.read_csv(input_dir / meta_df_csv)

    cov = _clinical_covariates(meta_df, index_col="File_ID")

    wide = meta_filtered.pivot_table(
        index="Sample", columns="Species", values="relative_abundance",
        fill_value=0)
    wide.index = wide.index.str.replace("_profile", "", regex=False)

    common = wide.index.intersection(cov.index)
    abundance = wide.loc[common]
    group = cov.loc[common, "sarc_status_bin"].map({0: "NonSarc", 1: "Sarc"})
    return abundance, group


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--meta-filtered-csv", required=True)
    ap.add_argument("--meta-df-csv", required=True)
    ap.add_argument("--alpha-metric", choices=list(ALPHA_METRICS),
                     default="shannon")
    ap.add_argument("--n-permutations", type=int, default=999)
    ap.add_argument("--out", default="diversity_panels")
    args = ap.parse_args()

    abundance, group = load_abundance(args.input_dir, args.meta_filtered_csv,
                                      args.meta_df_csv)

    # --- Alpha diversity ---
    metric_fn = ALPHA_METRICS[args.alpha_metric]
    alpha = abundance.apply(lambda row: metric_fn(row.values), axis=1)
    alpha_df = pd.DataFrame({"alpha": alpha, "group": group})

    nonsarc_vals = alpha_df.loc[alpha_df["group"] == "NonSarc", "alpha"]
    sarc_vals = alpha_df.loc[alpha_df["group"] == "Sarc", "alpha"]
    _, alpha_p = mannwhitneyu(nonsarc_vals, sarc_vals, alternative="two-sided")

    # --- Beta diversity ---
    dist = squareform(pdist(abundance.values, metric="braycurtis"))
    coords, var_explained = pcoa(dist, n_axes=2)
    pcoa_df = pd.DataFrame({"PC1": coords[:, 0], "PC2": coords[:, 1],
                            "group": group.values}, index=abundance.index)

    f_obs, r2_obs, perm_p = permanova(dist, group.values,
                                      n_perm=args.n_permutations)

    # --- Figure ---
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 5))

    # Panel A: alpha diversity boxplot
    box_data = [alpha_df.loc[alpha_df["group"] == g, "alpha"].values
               for g in GROUP_ORDER]
    bp = ax_a.boxplot(box_data, positions=range(len(GROUP_ORDER)),
                      widths=0.5, patch_artist=True, showfliers=False,
                      medianprops=dict(color="black", linewidth=1.4),
                      boxprops=dict(edgecolor="black", linewidth=0.9),
                      whiskerprops=dict(color="black", linewidth=0.9),
                      capprops=dict(color="black", linewidth=0.9))
    for patch, g in zip(bp["boxes"], GROUP_ORDER):
        patch.set_facecolor(GROUP_COLORS[g])
        patch.set_alpha(0.85)

    rng = np.random.default_rng(0)
    for i, g in enumerate(GROUP_ORDER):
        vals = alpha_df.loc[alpha_df["group"] == g, "alpha"].values
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax_a.scatter(np.full(len(vals), i) + jitter, vals, color="black",
                    s=16, alpha=0.7, zorder=3, edgecolor="white",
                    linewidth=0.4)

    ax_a.set_xticks(range(len(GROUP_ORDER)))
    ax_a.set_xticklabels(GROUP_ORDER, fontsize=11, fontweight="bold")
    ax_a.set_ylabel(f"{args.alpha_metric.capitalize()} diversity", fontsize=11)
    ax_a.set_title("A. Alpha diversity", fontsize=13, fontweight="bold", loc="left")
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    p_label = "p < 0.001" if alpha_p < 0.001 else f"p = {alpha_p:.3f}"
    y_max = alpha_df["alpha"].max()
    y_span = alpha_df["alpha"].max() - alpha_df["alpha"].min()
    bar_y = y_max + 0.06 * y_span
    ax_a.plot([0, 0, 1, 1], [bar_y, bar_y + 0.02 * y_span,
                              bar_y + 0.02 * y_span, bar_y],
             color="black", linewidth=1)
    ax_a.text(0.5, bar_y + 0.03 * y_span, p_label, ha="center", fontsize=10)
    ax_a.set_ylim(top=bar_y + 0.12 * y_span)

    # Panel B: PCoA
    for g in GROUP_ORDER:
        sub = pcoa_df[pcoa_df["group"] == g]
        ax_b.scatter(sub["PC1"], sub["PC2"], color=GROUP_COLORS[g],
                    label=g, s=32, alpha=0.85, edgecolor="black", linewidth=0.4)
        confidence_ellipse(ax_b, sub["PC1"].values, sub["PC2"].values,
                          GROUP_COLORS[g])

    ax_b.axhline(0, color="grey", linewidth=0.6, linestyle="--")
    ax_b.axvline(0, color="grey", linewidth=0.6, linestyle="--")
    ax_b.set_xlabel(f"PCo1 ({var_explained[0]*100:.1f}%)", fontsize=11)
    ax_b.set_ylabel(f"PCo2 ({var_explained[1]*100:.1f}%)", fontsize=11)
    ax_b.set_title("B. Beta diversity (Bray-Curtis PCoA)", fontsize=13,
                  fontweight="bold", loc="left")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.legend(frameon=False, fontsize=10, loc="best")

    perm_p_label = "p < 0.001" if perm_p < 0.001 else f"p = {perm_p:.3f}"
    ax_b.text(0.02, 0.98,
             f"PERMANOVA: R\u00b2 = {r2_obs:.3f}, {perm_p_label}",
             transform=ax_b.transAxes, fontsize=9, va="top",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))

    plt.tight_layout()
    fig.savefig(f"{args.out}.pdf", bbox_inches="tight")
    fig.savefig(f"{args.out}.png", dpi=300, bbox_inches="tight")
    print(f"Saved {args.out}.pdf and {args.out}.png")
    print(f"Alpha diversity ({args.alpha_metric}) Mann-Whitney p = {alpha_p:.4g}")
    print(f"PERMANOVA: pseudo-F = {f_obs:.3f}, R2 = {r2_obs:.3f}, p = {perm_p:.4g} "
         f"({args.n_permutations} permutations)")


if __name__ == "__main__":
    main()