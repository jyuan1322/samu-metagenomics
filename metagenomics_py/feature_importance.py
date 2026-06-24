"""feature_importance.py

Post-hoc feature-importance extraction from saved nested-CV results, shared
across assays. Each function reads a `cv_results` pickle (produced by
ml_engine.run_nested_cv) and summarizes a model's per-fold importances:

  * lasso_coef_summary  - L1 logistic-regression coefficients across folds
  * boruta_rf_summary   - Boruta selection frequency / ranking on the RF folds
  * fcnn_ig_summary     - FC_NN Integrated-Gradients attributions across folds

A helper, strip_group_prefix, removes the ColumnTransformer group prefix
(e.g. "clr_taxa__") that sklearn prepends to feature names.
"""
import numpy as np
import pandas as pd


def _fold_feature_names(cv_results):
    """Feature names out of the first fold's fitted ColumnTransformer."""
    return (cv_results["estimator"][0]
            .best_estimator_.named_steps["scaler"].get_feature_names_out())


def strip_group_prefix(series, prefixes=("clr_taxa__", "cont__", "bin__",
                                         "nmr__", "gcms__", "proteins__")):
    """Remove a leading ColumnTransformer group prefix from feature names."""
    out = series.astype(str)
    for p in prefixes:
        out = out.str.replace(p, "", regex=False)
    return out


# ---------------------------------------------------------------------------
# L1 logistic regression
# ---------------------------------------------------------------------------
def lasso_coef_summary(cv_results, outfile=None):
    """Summarize L1 logistic-regression coefficients across outer folds."""
    coefs = np.vstack([
        est.best_estimator_.named_steps["classifier"].coef_.ravel()
        for est in cv_results["estimator"]
    ])
    features = _fold_feature_names(cv_results)

    summary = pd.DataFrame({
        "feature": features,
        "mean_coef": coefs.mean(axis=0),
        "std_coef": coefs.std(axis=0),
        "nonzero_fraction": (coefs != 0).mean(axis=0),
        "lower_perc12_5": np.percentile(coefs, 12.5, axis=0),
        "upper_perc87_5": np.percentile(coefs, 87.5, axis=0),
    })
    summary["odds_ratio"] = np.exp(summary["mean_coef"])
    summary = summary.sort_values("nonzero_fraction", ascending=False)
    if outfile is not None:
        summary.to_csv(outfile, index=False)
    return summary


# ---------------------------------------------------------------------------
# Random forest via Boruta
# ---------------------------------------------------------------------------
def boruta_rf_summary(cv_results, X, y, outfile=None,
                      max_iter=100, base_random_state=42):
    """Run Boruta on each fold's trained RF and summarize selection/ranking.

    Boruta is fit per fold on that fold's scaled training data, using the
    fold's tuned RF as the base estimator. Requires the `boruta` package.
    """
    # Boruta references removed numpy aliases; restore them before import.
    np.int = int
    np.float = float
    np.bool = bool
    from boruta import BorutaPy

    rankings, selected, selected_weak = [], [], []
    for fold_idx, est in enumerate(cv_results["estimator"]):
        best_model = est.best_estimator_
        rf = best_model.named_steps["classifier"]
        train_idx = cv_results["indices"]["train"][fold_idx]

        X_train = best_model.named_steps["scaler"].transform(X.iloc[train_idx])
        y_train = y.iloc[train_idx]

        boruta = BorutaPy(rf, n_estimators="auto", max_iter=max_iter,
                          random_state=base_random_state + fold_idx)
        boruta.fit(X_train, y_train)
        rankings.append(boruta.ranking_)
        selected.append(boruta.support_)
        selected_weak.append(boruta.support_weak_)

    rankings = np.vstack(rankings)
    selected = np.vstack(selected)
    selected_weak = np.vstack(selected_weak)

    summary = pd.DataFrame({
        "feature": _fold_feature_names(cv_results),
        "selection_frequency": selected.mean(axis=0),
        "weak_selection_frequency": selected_weak.mean(axis=0),
        "mean_ranking": rankings.mean(axis=0),
        "median_ranking": np.median(rankings, axis=0),
        "min_ranking": rankings.min(axis=0),
        "max_ranking": rankings.max(axis=0),
    }).sort_values(["selection_frequency", "mean_ranking"],
                   ascending=[False, True])
    if outfile is not None:
        summary.to_csv(outfile, index=False)
    return summary


# ---------------------------------------------------------------------------
# FC_NN via Integrated Gradients
# ---------------------------------------------------------------------------
def fcnn_ig_summary(cv_results, X, outfile=None):
    """Integrated-Gradients attributions for the FC_NN, across outer folds.

    Attributions are computed on each fold's held-out test samples, on the
    class-1 logit. Requires torch and captum.
    """
    import torch
    from captum.attr import IntegratedGradients

    all_importances = []
    for fold_idx, est in enumerate(cv_results["estimator"]):
        pipeline = est.best_estimator_
        test_idx = cv_results["indices"]["test"][fold_idx]
        X_test = X.iloc[test_idx].astype(np.float32)

        X_input = pipeline.named_steps["scaler"].transform(X_test)
        X_input = torch.tensor(X_input, dtype=torch.float32)

        torch_model = pipeline.named_steps["classifier"].module_
        torch_model.eval()

        ig = IntegratedGradients(lambda x: torch_model(x)[:, 1])  # class-1 logit
        attributions, _ = ig.attribute(X_input, target=None,
                                       return_convergence_delta=True)
        all_importances.append(attributions.detach().numpy())

    all_importances = np.vstack(all_importances)
    features = X.columns if hasattr(X, "columns") else _fold_feature_names(cv_results)

    summary = pd.DataFrame({
        "feature": _fold_feature_names(cv_results),
        "mean_importance": all_importances.mean(axis=0),
        "std_importance": all_importances.std(axis=0),
        "positive_fraction": (all_importances > 0).mean(axis=0),
        "lower_perc12_5": np.percentile(all_importances, 12.5, axis=0),
        "upper_perc87_5": np.percentile(all_importances, 87.5, axis=0),
    }).sort_values("mean_importance", ascending=False)
    if outfile is not None:
        summary.to_csv(outfile, index=False)
    return summary


# ---------------------------------------------------------------------------
# Merge model summaries + DESeq2 into one table
# ---------------------------------------------------------------------------
def merge_importance_tables(deseq2_res, linreg_res, rf_res, fcnn_res,
                            outfile=None):
    """Outer-merge DESeq2, lasso, Boruta-RF, and FCNN summaries on `feature`.

    Each input is a DataFrame; the taxa group prefix is stripped from the model
    tables so their feature names line up with the DESeq2 feature column.
    """
    deseq2 = deseq2_res.rename(columns={"Unnamed: 0": "feature"})
    for df in (linreg_res, rf_res, fcnn_res):
        df["feature"] = strip_group_prefix(df["feature"])

    merged = pd.merge(
        deseq2[["feature", "log2FoldChange", "padj"]],
        linreg_res[["feature", "mean_coef", "std_coef", "odds_ratio",
                    "nonzero_fraction", "lower_perc12_5", "upper_perc87_5"]],
        on="feature", how="outer",
    ).rename(columns={
        "log2FoldChange": "deseq2_log2FC", "padj": "deseq2_padj",
        "mean_coef": "linreg_mean_coef", "std_coef": "linreg_std_coef",
        "odds_ratio": "linreg_odds_ratio",
        "nonzero_fraction": "linreg_nonzero_fraction",
        "lower_perc12_5": "linreg_lower_perc12_5",
        "upper_perc87_5": "linreg_upper_perc87_5",
    })

    merged = pd.merge(
        merged,
        rf_res[["feature", "selection_frequency", "weak_selection_frequency",
                "mean_ranking", "median_ranking", "min_ranking", "max_ranking"]],
        on="feature", how="outer",
    ).rename(columns={
        "selection_frequency": "rf_selection_frequency",
        "weak_selection_frequency": "rf_weak_selection_frequency",
        "mean_ranking": "rf_mean_ranking", "median_ranking": "rf_median_ranking",
        "min_ranking": "rf_min_ranking", "max_ranking": "rf_max_ranking",
    })

    merged = pd.merge(
        merged,
        fcnn_res[["feature", "mean_importance", "std_importance",
                  "positive_fraction", "lower_perc12_5", "upper_perc87_5"]],
        on="feature", how="outer",
    ).rename(columns={
        "mean_importance": "fcnn_mean_importance",
        "std_importance": "fcnn_std_importance",
        "positive_fraction": "fcnn_positive_fraction",
        "lower_perc12_5": "fcnn_lower_perc12_5",
        "upper_perc87_5": "fcnn_upper_perc87_5",
    })

    if outfile is not None:
        merged.to_csv(outfile, index=False)
    return merged
