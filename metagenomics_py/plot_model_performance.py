"""
plot_model_performance.py

Publication-quality figure: outer-fold cross-validated performance, grouped
by model and by which feature set was used (all features / covariates only /
microbiome only), shown as grouped box-and-whisker plots with individual
fold scores overlaid as dots.

Reads the cv_results_<model>_<subset>.pkl files written by
ml_engine.run_model_set() (metagenomics_py), pulls cv_results["test_score"]
(one value per outer fold) out of each, and plots grouped bars (mean across
folds) with the individual fold scores overlaid as dots.

IMPORTANT: cv_results["test_score"] only reflects the scoring metric you
actually asked for (e.g. roc_auc) if ml_engine.py's cross_validate() call
passes scoring=scoring_metric. If your pickles predate that fix, this script
will be plotting accuracy, not AUC -- check --metric-label against what was
actually saved.

Expected filenames (as written by run_model_set):
  cv_results_<model>_full.pkl
  cv_results_<model>_<feature_group_name>_only.pkl   (e.g. clr_taxa_only)
  cv_results_<model>_clinical_only.pkl

Usage:
  python plot_model_performance.py <results_dir> \
      --models lasso_logreg random_forest fcnn \
      --feature-group clr_taxa \
      --metric-label "ROC AUC" \
      --out model_performance
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

MODEL_DISPLAY = {
    "lasso_logreg": "Lasso logistic\nregression",
    "random_forest": "Random forest",
    "fcnn": "FCNN",
    "mlp": "MLP",
    "adaboost": "AdaBoost",
}

SUBSET_DISPLAY = {
    "full": "All features",
    "clinical_only": "Covariates only",
    # the microbiome-only key is built from --feature-group at runtime,
    # e.g. "clr_taxa_only" -> "Microbiome only"
}

SUBSET_COLORS = {
    "All features": "#4C72B0",
    "Covariates only": "#DD8452",
    "Microbiome only": "#55A868",
}


def load_fold_scores(results_dir, models, feature_group, metric_key="test_score"):
    results_dir = Path(results_dir)
    subset_display = dict(SUBSET_DISPLAY)
    subset_display[f"{feature_group}_only"] = "Microbiome only"

    rows = []
    missing = []
    for model in models:
        for subset_key, subset_label in subset_display.items():
            path = results_dir / f"cv_results_{model}_{subset_key}.pkl"
            if not path.exists():
                missing.append(path.name)
                continue
            cv_results = pd.read_pickle(path)
            scores = np.asarray(cv_results[metric_key])
            for fold_i, score in enumerate(scores, 1):
                rows.append({
                    "model": MODEL_DISPLAY.get(model, model),
                    "subset": subset_label,
                    "fold": fold_i,
                    "score": score,
                })

    if missing:
        print(f"Warning: {len(missing)} expected file(s) not found, skipped: "
              f"{missing}")

    if not rows:
        raise FileNotFoundError(
            "No cv_results files found -- check results_dir/models/feature_group.")

    return pd.DataFrame(rows)


def plot_performance(df, metric_label="ROC AUC", out_prefix="model_performance"):
    models = list(dict.fromkeys(df["model"]))  # preserve first-seen order
    subsets = [s for s in SUBSET_COLORS if s in df["subset"].unique()]

    n_model = len(models)
    n_subset = len(subsets)
    box_w = 0.8 / n_subset
    x = np.arange(n_model)

    fig, ax = plt.subplots(figsize=(1.8 * n_model + 2, 5))

    rng = np.random.default_rng(0)
    for si, subset in enumerate(subsets):
        offset = (si - (n_subset - 1) / 2) * box_w
        positions, box_data = [], []
        for mi, model in enumerate(models):
            vals = df[(df["model"] == model) & (df["subset"] == subset)]["score"]
            positions.append(x[mi] + offset)
            box_data.append(vals.values)

            # individual fold points, jittered slightly within the box width
            jitter = rng.uniform(-box_w * 0.25, box_w * 0.25, size=len(vals))
            ax.scatter(np.full(len(vals), x[mi] + offset) + jitter, vals,
                       color="black", s=14, zorder=3, alpha=0.7,
                       edgecolor="white", linewidth=0.4)

        bp = ax.boxplot(box_data, positions=positions, widths=box_w * 0.8,
                        patch_artist=True, showfliers=False, zorder=2,
                        medianprops=dict(color="black", linewidth=1.2),
                        boxprops=dict(edgecolor="black", linewidth=0.8),
                        whiskerprops=dict(color="black", linewidth=0.8),
                        capprops=dict(color="black", linewidth=0.8))
        for patch in bp["boxes"]:
            patch.set_facecolor(SUBSET_COLORS[subset])
            patch.set_alpha(0.85)

    legend_handles = [Patch(facecolor=SUBSET_COLORS[s], edgecolor="black",
                            alpha=0.85, label=s) for s in subsets]

    ax.set_xticks(x)
    ax.set_xticklabels(models, fontsize=10, fontweight="bold")
    ax.set_ylabel(metric_label, fontsize=11)
    ax.set_title(f"Model performance by feature set ({metric_label}, outer CV folds)",
                fontsize=12, fontweight="bold")
    ax.legend(handles=legend_handles, title="Feature set", loc="upper center",
              bbox_to_anchor=(0.5, -0.12), ncol=n_subset, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.6, n_model - 0.4)
    plt.tight_layout()

    fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches="tight")
    print(f"Saved {out_prefix}.pdf and {out_prefix}.png")
    return fig


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", nargs="?", default=None,
                     help="Directory containing cv_results_*.pkl "
                          "(not needed if --from-csv is used)")
    ap.add_argument("--models", nargs="+",
                     default=["lasso_logreg", "random_forest", "fcnn"])
    ap.add_argument("--feature-group", default="clr_taxa",
                     help="feature_group_name used when training "
                          "(the *_only.pkl prefix), e.g. clr_taxa")
    ap.add_argument("--metric-key", default="test_score",
                     help="Key inside the cv_results dict to plot")
    ap.add_argument("--metric-label", default="ROC AUC",
                     help="Axis/title label -- must match what --metric-key "
                          "actually contains")
    ap.add_argument("--from-csv", default=None,
                     help="Path to a long-format CSV with columns "
                          "model,subset,fold,score (e.g. from "
                          "extract_fold_auc.py) -- bypasses reading "
                          "cv_results pickles directly.")
    ap.add_argument("--out", default="model_performance")
    args = ap.parse_args()

    if args.from_csv:
        df = pd.read_csv(args.from_csv)
    else:
        df = load_fold_scores(args.results_dir, args.models, args.feature_group,
                              metric_key=args.metric_key)
    plot_performance(df, metric_label=args.metric_label, out_prefix=args.out)