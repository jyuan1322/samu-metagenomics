"""ml_engine.py

Shared nested cross-validation engine for the SaMu sarcopenia-prediction
models. Every assay entry script (metagenomics, mass-spec, NMR, LC/GC-MS)
imports `run_nested_cv` from here; only the data loading differs between them.

The engine builds a preprocessing ColumnTransformer from a set of named feature
groups, wraps it with a classifier in a Pipeline, and runs nested CV: an inner
GridSearchCV for hyperparameter tuning inside an outer StratifiedKFold for
unbiased performance estimation. It saves the per-fold `cross_validate` results
(estimators + indices, for downstream feature-importance extraction) and a
mean-ROC plot.
"""
import pickle
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import (
    StratifiedKFold, GridSearchCV, cross_val_score, cross_validate,
)
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc, RocCurveDisplay
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier

from fcnn_wrapper import make_ffnn_classifier


# ---------------------------------------------------------------------------
# Hyperparameter grids
# ---------------------------------------------------------------------------
def _lasso_C_path(lambda_min=0.001, l_max=100, path_len=100):
    """Regularization-strength path for L1 logistic regression.

    Mirrors the lambda path used in Dawkins et al. 2024 (mmethane benchmarker):
    a log-spaced lambda grid converted to sklearn's C = 1 / lambda.
    """
    l_path = np.logspace(np.log10(l_max * lambda_min), np.log10(l_max), path_len)
    return [1.0 / l for l in l_path]


def _ensemble_param_grid(model_type):
    """Build a classifier__-prefixed grid for RF / AdaBoost, keeping only the
    keys each estimator actually accepts."""
    grid = {
        "bootstrap": [True],
        "n_estimators": [50, 100],
        "max_depth": [None],
        "max_features": [None, "sqrt"],
        "min_samples_split": [2, 9],
        "min_samples_leaf": [1, 5],
        "learning_rate": [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3,
                          1e-2, 5e-2, 0.1, 0.5, 1, 5, 10],
    }
    estimator = (RandomForestClassifier() if model_type == "random_forest"
                 else AdaBoostClassifier())
    allowed = set(estimator._parameter_constraints.keys())
    grid = {k: v for k, v in grid.items() if k in allowed}
    return {f"classifier__{k}": v for k, v in grid.items()}


def build_model(model_type, fcnn_input_dim=None):
    """Return (estimator, param_grid) for a given model_type.

    model_type is one of: lasso_logreg, mlp, fcnn, random_forest, adaboost.
    fcnn_input_dim is required for the fcnn estimator (set after the
    ColumnTransformer output width is known).
    """
    if model_type in ("random_forest", "adaboost"):
        param_grid = _ensemble_param_grid(model_type)
    else:
        param_grid = None

    model_dict = {
        "lasso_logreg": (
            LogisticRegression(penalty="l1", solver="saga",
                               max_iter=5000, random_state=42),
            {"classifier__C": _lasso_C_path()},
        ),
        "mlp": (
            MLPClassifier(max_iter=1000, random_state=42),
            {"classifier__hidden_layer_sizes": [(50,), (100,), (50, 50)],
             "classifier__alpha": [1e-4, 1e-3, 1e-2]},
        ),
        "fcnn": (
            make_ffnn_classifier(input_dim=fcnn_input_dim, h=[0, 0]),
            {"classifier__optimizer__weight_decay": [0.01],
             "classifier__module__p": [0.1],
             "classifier__lr": [0.001],
             "classifier__max_epochs": [2000]},
        ),
        "random_forest": (RandomForestClassifier(random_state=42), param_grid),
        "adaboost": (AdaBoostClassifier(random_state=42), param_grid),
    }
    if model_type not in model_dict:
        raise ValueError(f"Unknown model_type '{model_type}'. "
                         f"Choose from {list(model_dict)}")
    return model_dict[model_type]


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------
def build_column_transformer(feature_groups, binary_cols=None):
    """Build a ColumnTransformer from named feature groups.

    feature_groups : dict {name: [columns]} for groups that should be
        StandardScaler-scaled (e.g. {"cont": [...], "clr_taxa": [...]}).
        Empty/None groups are skipped.
    binary_cols : columns passed through unscaled.

    Returns the ColumnTransformer, or the string 'passthrough' if no group
    applies (so the pipeline still runs on a covariate-free matrix).
    """
    transformers = []
    for name, cols in (feature_groups or {}).items():
        if cols:
            transformers.append((name, StandardScaler(), list(cols)))
    if binary_cols:
        transformers.append(("bin", "passthrough", list(binary_cols)))
    return ColumnTransformer(transformers=transformers) if transformers else "passthrough"


# ---------------------------------------------------------------------------
# Nested cross-validation
# ---------------------------------------------------------------------------
def run_nested_cv(
        X, y,
        model_type="lasso_logreg",
        feature_groups=None,
        binary_cols=None,
        outfile_cv="cv_results.pkl",
        outfile_roc_curve="roc_curve.pdf",
        outer_seed=222,
        inner_seed=111,
        scoring_metric="roc_auc"):
    """Run nested CV for one model on one feature matrix.

    feature_groups : dict {name: [columns]} of scaled feature groups (taxa,
        metabolites, proteins, continuous covariates, ...). binary_cols are
        passed through unscaled. Together they define the ColumnTransformer.

    Saves the cross_validate result (with estimators + fold indices) to
    outfile_cv and a mean-ROC plot to outfile_roc_curve.
    """
    col_scaler = build_column_transformer(feature_groups, binary_cols)

    # FC_NN needs its input dimension = transformed feature width.
    if col_scaler == "passthrough":
        fcnn_input_dim = X.shape[1]
    else:
        fcnn_input_dim = col_scaler.fit_transform(X).shape[1]
    print(f"Input dimension after preprocessing: {fcnn_input_dim}")

    estimator, param_grid = build_model(model_type, fcnn_input_dim=fcnn_input_dim)
    pipeline = Pipeline([("scaler", col_scaler), ("classifier", estimator)])

    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=outer_seed)
    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=inner_seed)

    inner_grid_search = GridSearchCV(
        estimator=pipeline, param_grid=param_grid,
        cv=inner_cv, scoring=scoring_metric, n_jobs=-1,
    )

    nested_scores = cross_val_score(
        inner_grid_search, X=X, y=y, cv=outer_cv,
        scoring=scoring_metric, n_jobs=-1,
    )
    print(f"{scoring_metric} (outer folds): {nested_scores}")

    cv_results = cross_validate(
        inner_grid_search, X=X, y=y, cv=outer_cv,
        return_estimator=True, return_indices=True, n_jobs=-1,
    )
    with open(outfile_cv, "wb") as f:
        pickle.dump(cv_results, f)

    _plot_mean_roc(cv_results, X, y, outer_cv, outfile_roc_curve)
    return cv_results


def run_model_set(
        X, y, feature_groups, binary_cols, output_dir,
        models=("lasso_logreg", "random_forest", "fcnn"),
        feature_group_name=None,
        outer_seed=222, inner_seed=111, scoring_metric="roc_auc"):
    """Run the standard full / feature-only / clinical-only sweep for each model.

    For each model in `models`, three nested-CV runs are executed:
      * full           - covariates + the assay feature group(s)
      * <feature>_only - assay feature group(s) only (no covariates)
      * clinical_only  - covariates only

    feature_group_name names the assay feature group used for the *_only run
    (e.g. "clr_taxa", "proteins", "gcms"); if None it is inferred as the single
    non-continuous, non-binary group in feature_groups.

    Output files are written into output_dir as cv_results_<model>_<subset>.pkl
    and roc_curve_<model>_<subset>.pdf. FC_NN runs on a float32 copy of X.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cont = {k: v for k, v in feature_groups.items() if k == "cont"}
    assay = {k: v for k, v in feature_groups.items() if k != "cont"}
    if feature_group_name is None and assay:
        feature_group_name = next(iter(assay))

    subsets = {
        "full": (feature_groups, binary_cols),
        f"{feature_group_name}_only": (assay, None),
        "clinical_only": ({**cont}, binary_cols),
    }

    for model in models:
        X_model = X.astype(np.float32) if model == "fcnn" else X
        for subset, (groups, bins) in subsets.items():
            run_nested_cv(
                X_model, y, model_type=model,
                feature_groups=groups, binary_cols=bins,
                outfile_cv=output_dir / f"cv_results_{model}_{subset}.pkl",
                outfile_roc_curve=output_dir / f"roc_curve_{model}_{subset}.pdf",
                outer_seed=outer_seed, inner_seed=inner_seed,
                scoring_metric=scoring_metric,
            )


def _plot_mean_roc(cv_results, X, y, outer_cv, outfile_roc_curve):
    """Plot per-fold ROC curves plus the interpolated mean +/- 1 s.d. band.

    Follows the scikit-learn ROC-with-CV example.
    """
    n_splits = outer_cv.get_n_splits()
    mean_fpr = np.linspace(0, 1, 100)
    interp_tprs = []

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    curve_kwargs_list = [dict(alpha=0.3, lw=1, color=colors[i % len(colors)])
                         for i in range(n_splits)]
    names = [f"ROC fold {i + 1}" for i in range(n_splits)]

    _, ax = plt.subplots(figsize=(6, 6))
    viz = RocCurveDisplay.from_cv_results(
        cv_results, X, y, ax=ax, name=names,
        curve_kwargs=curve_kwargs_list, plot_chance_level=True,
    )

    for idx in range(n_splits):
        interp_tpr = np.interp(mean_fpr, viz.fpr[idx], viz.tpr[idx])
        interp_tpr[0] = 0.0
        interp_tprs.append(interp_tpr)

    mean_tpr = np.mean(interp_tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(viz.roc_auc)

    ax.plot(mean_fpr, mean_tpr, color="b", lw=2, alpha=0.8,
            label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc))

    std_tpr = np.std(interp_tprs, axis=0)
    ax.fill_between(
        mean_fpr,
        np.maximum(mean_tpr - std_tpr, 0),
        np.minimum(mean_tpr + std_tpr, 1),
        color="grey", alpha=0.2, label=r"$\pm$ 1 std. dev.",
    )

    ax.set(xlabel="False Positive Rate", ylabel="True Positive Rate",
           title="Mean ROC curve")
    ax.legend(loc="lower right")
    plt.savefig(outfile_roc_curve, format="pdf", bbox_inches="tight")
    plt.close()
