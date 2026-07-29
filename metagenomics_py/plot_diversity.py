"""
plot_diversity.py

Two-panel publication figure:
  (A) Alpha diversity (Shannon index by default) per sample, as a boxplot
      split by Sarc / NoSarc, with individual points overlaid and a
      Mann-Whitney U p-value annotated. Computed independently in Python
      from raw MetaPhlAn relative abundances (same filtering as the R
      pipeline, via loaders._clinical_covariates).
  (B) Beta diversity: PCoA (Bray-Curtis) scatter, colored by group. This
      panel does NOT recompute PCoA in Python -- it reads the exact
      coordinates and % variance explained that 02_diversity.R exported
      right after its own ordinate(ps, method="PCoA", ...) call
      (beta_diversity_pcoa_coords.csv, beta_diversity_pcoa_variance_explained.csv),
      so Panel B is guaranteed to be pixel-for-pixel the same ordination as
      R's beta_diversity_group.pdf rather than a second, independently
      re-derived one (avoiding any mismatch in eigen-decomposition
      convention, e.g. how ape::pcoa's negative eigenvalues from this
      non-Euclidean Bray-Curtis distance get folded into % variance
      explained).

Uses the SAME sample filtering as the rest of the pipeline (Full.SaMu==1,
age_def>=50, complete covariates) via loaders._clinical_covariates for
Panel A. Panel B is driven entirely by the R-exported CSVs -- it does not
independently load or filter the abundance table.

Usage:
  python plot_diversity.py \
      --input-dir /data/local/jy1008/SaMu/results/latest/metagenomics_R \
      --meta-filtered-csv meta_filtered_06242026.csv \
      --meta-df-csv meta_df_FullSaMu_06242026.csv \
      --pcoa-coords-csv beta_diversity_pcoa_coords.csv \
      --pcoa-variance-csv beta_diversity_pcoa_variance_explained.csv \
      --alpha-metric shannon \
      --n-permutations 999 \
      --out diversity_panels
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu

from loaders import _clinical_covariates

# GROUP_COLORS = {"NoSarc": "#4C72B0", "Sarc": "#C44E52"}
GROUP_COLORS = {"NoSarc": "#4DAF4A", "Sarc": "#984EA3"}
GROUP_ORDER = ["NoSarc", "Sarc"]

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]  # falls back to DejaVu Sans if Arial isn't installed

# ---------------------------------------------------------------------------
# Alpha diversity metrics (computed per sample from a relative-abundance row)
# ---------------------------------------------------------------------------
def shannon(p):
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())

def shannon_normalized(row_values, log_total_n_taxa):
    p = row_values / row_values.sum()      # normalize to proportions, like vegan::diversity
    p = p[p > 0]
    h = -(p * np.log(p)).sum()
    return h / log_total_n_taxa             # divide by log(total taxa in table), not per-sample richness

def simpson(p):
    return float(1 - (p ** 2).sum())


def richness(p):
    return float((p > 0).sum())


ALPHA_METRICS = {"normalized shannon": shannon_normalized, "simpson": simpson, "richness": richness}


# ---------------------------------------------------------------------------
# PERMANOVA (Anderson 2001 pseudo-F, permutation p-value). Still computed
# independently in Python (printed/returned, matching adonis2's console
# output in 02_diversity.R) -- only the PCoA *coordinates* are imported from
# R now, not the PERMANOVA test itself.
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
    group = cov.loc[common, "sarc_status_bin"].map({0: "NoSarc", 1: "Sarc"})
    return abundance, group


def load_pcoa_from_r(input_dir, coords_csv, variance_csv, group_col):
    """Read the PCoA coordinates + % variance explained that
    02_diversity.R exported right after its ordinate() call, so Panel B
    plots the SAME numbers R did rather than an independently re-derived
    ordination. coords_csv is expected to have columns File_ID, PC1, PC2,
    and <group_col> (whatever GROUP_VAR was in the R config); group values
    are mapped 0/1 -> NoSarc/Sarc the same way Panel A's group labels are,
    unless the R column already contains those string labels."""
    input_dir = Path(input_dir)
    coords = pd.read_csv(input_dir / coords_csv)
    variance = pd.read_csv(input_dir / variance_csv).set_index("axis")["var_explained"]

    # A single canonical label map, used whether the R column is numeric
    # (0/1) or already text -- if R exported a differently-spelled label
    # (e.g. "NonSarc") this normalizes it to match Panel A instead of
    # passing it through as-is, which is what silently caused the two
    # panels to disagree before.
    label_map = {0: "NoSarc", 1: "Sarc", "NonSarc": "NoSarc", "NoSarc": "NoSarc",
                "Sarc": "Sarc"}
    coords["group"] = coords[group_col].map(label_map).fillna(coords[group_col])

    var_explained = np.array([variance.loc["PC1"], variance.loc["PC2"]])
    return coords.set_index("File_ID"), var_explained


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--meta-filtered-csv", required=True)
    ap.add_argument("--meta-df-csv", required=True)
    ap.add_argument("--pcoa-coords-csv", default="beta_diversity_pcoa_coords.csv",
                     help="PCoA coordinates exported by 02_diversity.R right "
                          "after its ordinate() call (default: "
                          "beta_diversity_pcoa_coords.csv)")
    ap.add_argument("--pcoa-variance-csv",
                     default="beta_diversity_pcoa_variance_explained.csv",
                     help="% variance explained per axis, exported by "
                          "02_diversity.R (default: "
                          "beta_diversity_pcoa_variance_explained.csv)")
    ap.add_argument("--pcoa-group-col", default="sarc_status_bin",
                     help="Group column name inside --pcoa-coords-csv (i.e. "
                          "whatever GROUP_VAR was set to in R's config.R; "
                          "default sarc_status_bin)")
    ap.add_argument("--alpha-metric", choices=list(ALPHA_METRICS),
                     default="Normalized shannon")
    ap.add_argument("--n-permutations", type=int, default=999)
    ap.add_argument("--out", default="diversity_panels")
    args = ap.parse_args()

    abundance, group = load_abundance(args.input_dir, args.meta_filtered_csv,
                                      args.meta_df_csv)

    # --- Alpha diversity (Panel A) ---
    # metric_fn = ALPHA_METRICS[args.alpha_metric]
    # alpha = abundance.apply(lambda row: metric_fn(row.values), axis=1)
    log_total_n_taxa = np.log(abundance.shape[1])   # == log(nrow(tax_table(ps))) in R
    alpha = abundance.apply(
        lambda row: shannon_normalized(row.values, log_total_n_taxa), axis=1)


    alpha_df = pd.DataFrame({"alpha": alpha, "group": group})

    nosarc_vals = alpha_df.loc[alpha_df["group"] == "NoSarc", "alpha"]
    sarc_vals = alpha_df.loc[alpha_df["group"] == "Sarc", "alpha"]
    _, alpha_p = mannwhitneyu(nosarc_vals, sarc_vals, alternative="two-sided")

    # --- Beta diversity (Panel B) -- read R's own PCoA output directly ---
    pcoa_df, var_explained = load_pcoa_from_r(
        args.input_dir, args.pcoa_coords_csv, args.pcoa_variance_csv,
        args.pcoa_group_col)

    # PERMANOVA is still computed in Python (for the console summary below)
    # from a fresh Bray-Curtis distance on the same abundance table -- only
    # the PCoA coordinates/variance-explained are taken from R now.
    common = abundance.index.intersection(pcoa_df.index)
    dist = squareform(pdist(abundance.loc[common].values, metric="braycurtis"))
    f_obs, r2_obs, perm_p = permanova(dist, group.loc[common].values,
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
    ax_a.set_title("Alpha diversity", fontsize=13, fontweight="bold", loc="left")
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

    # Panel B: PCoA -- plotted directly from R's exported coordinates.
    # Matches R's plot_ordination(..., color = GROUP_VAR) + geom_point(size = 3)
    # + theme_minimal(): points colored by group only, no ellipses, no
    # crosshairs. PERMANOVA R2/p ARE drawn here (unlike R's plot, which
    # leaves that to adonis_group_results.csv) -- added back per request.
    for g in GROUP_ORDER:
        sub = pcoa_df[pcoa_df["group"] == g]
        ax_b.scatter(sub["PC1"], sub["PC2"], color=GROUP_COLORS[g],
                    label=g, s=32, alpha=0.85, edgecolor="black", linewidth=0.4)

    ax_b.set_xlabel(f"PCo1 ({var_explained[0]*100:.1f}%)", fontsize=11)
    ax_b.set_ylabel(f"PCo2 ({var_explained[1]*100:.1f}%)", fontsize=11)
    ax_b.set_title("Beta diversity (Bray-Curtis PCoA)", fontsize=13,
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
         f"({args.n_permutations} permutations) -- also drawn on Panel B")


if __name__ == "__main__":
    main()