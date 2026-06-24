import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_score, cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay
from fcnn_wrapper import *

# ----------------------------------------
# Paths
# ----------------------------------------
experiment_name = "GC_MS"
input_dir = Path("/data/local/jy1008/SaMu/results/latest/proteomics_GC-MS")
output_dir = Path("/data/local/jy1008/SaMu/results/latest/proteomics_GC-MS/ml_models")
output_dir.mkdir(parents=True, exist_ok=True)

# ----------------------------------------
# Load GC-MS data (output from your R DEP2 pipeline)
# ----------------------------------------
gcms_df = pd.read_csv(input_dir / f"proteomics_{experiment_name}_vsn_imputed_matrix.csv")
meta_df = pd.read_csv(input_dir / f"proteomics_{experiment_name}_vsn_imputed_metadata.csv")

gcms_df = gcms_df.set_index("Sample")
gcms_df.index = "SaMu" + gcms_df.index.str.extract(r"(\d+)$")[0]
gcms_df.index.name = "Sample"
meta_df = meta_df.set_index("label")

# Drop unknown/missing diagnosis
meta_df = meta_df[meta_df["sarc_status_bin"].notna()]
meta_df = meta_df[meta_df["sarc_status_bin"] != "Unknown"]

# Binarize: Sarc=1, NoSarc=0
meta_df["sarc_status_bin"] = (meta_df["sarc_status_bin"] == "Sarc").astype(int)

# Binarize sex
meta_df["sex"] = meta_df["sex"].map({"M": 0, "F": 1})  # or whichever encoding you prefer

# ----------------------------------------
# Align samples between features and metadata
# ----------------------------------------
common_samples = gcms_df.index.intersection(meta_df.index)
print(f"Common samples: {len(common_samples)}")

gcms_df = gcms_df.loc[common_samples]
meta_df = meta_df.loc[common_samples]

# Drop samples with missing covariates
valid_samples = meta_df[["age_def", "sex", "smke", "alco", "nutr_score", "bmi"]].dropna().index
common_samples = common_samples.intersection(valid_samples)
print(f"Samples after dropping covariate NAs: {len(common_samples)}")

gcms_df = gcms_df.loc[common_samples]
meta_df = meta_df.loc[common_samples]

y = meta_df["sarc_status_bin"]

continuous_cols = ["age_def", "nutr_score", "bmi"]
binary_cols     = ["sex", "smke", "alco"]
gcms_cols       = gcms_df.columns.tolist()

# replace brackets and other characters in column names
gcms_df.columns = gcms_df.columns.str.replace(r"[^\w]", "_", regex=True).str.strip()
gcms_cols = gcms_df.columns.tolist()

X = pd.concat([meta_df[continuous_cols + binary_cols], gcms_df], axis=1)

assert "sarc_status_bin" not in X.columns
assert X.shape[0] == y.shape[0]
print(f"Feature matrix: {X.shape}, class balance:\n{y.value_counts()}")


# ----------------------------------------
# ML pipeline (same as NMR script)
# ----------------------------------------
def run_ml_pipeline(
        X, y,
        model_type="lasso_logreg",
        continuous_cols=None,
        binary_cols=None,
        gcms_cols=None,
        outfile_cv="cv_results.pkl",
        outfile_roc_curve="roc_curve.pdf"):

    transformers = []
    if continuous_cols:
        transformers.append(("cont", StandardScaler(), continuous_cols))
    if binary_cols:
        transformers.append(("bin", "passthrough", binary_cols))
    if gcms_cols:
        transformers.append(("gcms", StandardScaler(), gcms_cols))

    if transformers:
        col_scaler = ColumnTransformer(transformers=transformers)
        X_transformed = col_scaler.fit_transform(X)
    else:
        col_scaler = "passthrough"
        X_transformed = X.values

    fcnn_input_dim = X_transformed.shape[1]
    print(f"Input dim: {fcnn_input_dim}")

    # Lasso C path (same as NMR script)
    lambda_min = 0.001
    path_len   = 100
    l_max      = 100
    l_path     = np.logspace(np.log10(l_max * lambda_min), np.log10(l_max), path_len)
    Cs         = [1 / l for l in l_path]

    # Random forest param grid
    rf_param_grid = {
        "classifier__n_estimators":    [50, 100],
        "classifier__max_depth":       [None],
        "classifier__max_features":    [None, "sqrt"],
        "classifier__min_samples_split": [2, 9],
        "classifier__min_samples_leaf":  [1, 5],
    }

    model_dict = {
        "lasso_logreg": {
            "estimator":  LogisticRegression(penalty="l1", solver="saga",
                                             max_iter=5000, random_state=42),
            "param_grid": {"classifier__C": Cs},
        },
        "random_forest": {
            "estimator":  RandomForestClassifier(random_state=42),
            "param_grid": rf_param_grid,
        },
        "fcnn": {
            "estimator":  make_ffnn_classifier(input_dim=fcnn_input_dim, h=[0, 0]),
            "param_grid": {
                "classifier__optimizer__weight_decay": [0.01],
                "classifier__module__p":   [0.1],
                "classifier__lr":          [0.001],
                "classifier__max_epochs":  [2000],
            },
        },
    }

    if model_type not in model_dict:
        raise ValueError(f"Unknown model_type '{model_type}'")

    pipeline = Pipeline([
        ("scaler",     col_scaler),
        ("classifier", model_dict[model_type]["estimator"]),
    ])
    param_grid = model_dict[model_type]["param_grid"]

    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=222)
    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=111)

    inner_grid_search = GridSearchCV(
        estimator=pipeline,
        param_grid=param_grid,
        cv=inner_cv,
        scoring="roc_auc",
        n_jobs=-1,
    )

    nested_scores = cross_val_score(
        inner_grid_search, X=X, y=y,
        cv=outer_cv, scoring="roc_auc", n_jobs=-1,
    )
    print(f"Nested CV AUC: {nested_scores}  |  mean={nested_scores.mean():.3f}")

    cv_results = cross_validate(
        inner_grid_search, X=X, y=y,
        cv=outer_cv,
        return_estimator=True,
        return_indices=True,
        n_jobs=-1,
    )
    with open(outfile_cv, "wb") as f:
        pickle.dump(cv_results, f)

    # ROC plot
    n_splits = outer_cv.get_n_splits()
    mean_fpr   = np.linspace(0, 1, 100)
    interp_tprs = []

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    curve_kwargs_list = [dict(alpha=0.3, lw=1, color=colors[i % len(colors)])
                         for i in range(n_splits)]
    names = [f"ROC fold {i+1}" for i in range(n_splits)]

    _, ax = plt.subplots(figsize=(6, 6))
    viz = RocCurveDisplay.from_cv_results(
        cv_results, X, y, ax=ax,
        name=names, curve_kwargs=curve_kwargs_list, plot_chance_level=True,
    )
    for idx in range(n_splits):
        interp_tpr = np.interp(mean_fpr, viz.fpr[idx], viz.tpr[idx])
        interp_tpr[0] = 0.0
        interp_tprs.append(interp_tpr)

    mean_tpr      = np.mean(interp_tprs, axis=0);  mean_tpr[-1] = 1.0
    mean_auc      = auc(mean_fpr, mean_tpr)
    std_auc       = np.std(viz.roc_auc)
    std_tpr       = np.std(interp_tprs, axis=0)

    ax.plot(mean_fpr, mean_tpr, color="b", lw=2, alpha=0.8,
            label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc))
    ax.fill_between(mean_fpr,
                    np.maximum(mean_tpr - std_tpr, 0),
                    np.minimum(mean_tpr + std_tpr, 1),
                    color="grey", alpha=0.2, label=r"$\pm$ 1 std. dev.")
    ax.set(xlabel="False Positive Rate", ylabel="True Positive Rate",
           title=f"Mean ROC — {model_type} — {experiment_name}")
    ax.legend(loc="lower right")
    plt.savefig(outfile_roc_curve, format="pdf", bbox_inches="tight")
    plt.close()

# ----------------------------------------
# Run models
# ----------------------------------------

# Lasso logistic regression
run_ml_pipeline(X=X, y=y, model_type="lasso_logreg",
                continuous_cols=continuous_cols, binary_cols=binary_cols, gcms_cols=gcms_cols,
                outfile_cv=str(output_dir / "cv_results_l1logreg_gcms_and_covars.pkl"),
                outfile_roc_curve=str(output_dir / "roc_l1logreg_gcms_and_covars.pdf"))

run_ml_pipeline(X=X[gcms_cols], y=y, model_type="lasso_logreg",
                gcms_cols=gcms_cols,
                outfile_cv=str(output_dir / "cv_results_l1logreg_gcms_only.pkl"),
                outfile_roc_curve=str(output_dir / "roc_l1logreg_gcms_only.pdf"))

run_ml_pipeline(X=X[continuous_cols + binary_cols], y=y, model_type="lasso_logreg",
                continuous_cols=continuous_cols, binary_cols=binary_cols,
                outfile_cv=str(output_dir / "cv_results_l1logreg_clinical_only.pkl"),
                outfile_roc_curve=str(output_dir / "roc_l1logreg_clinical_only.pdf"))

# Random forest
run_ml_pipeline(X=X, y=y, model_type="random_forest",
                continuous_cols=continuous_cols, binary_cols=binary_cols, gcms_cols=gcms_cols,
                outfile_cv=str(output_dir / "cv_results_rf_gcms_and_covars.pkl"),
                outfile_roc_curve=str(output_dir / "roc_rf_gcms_and_covars.pdf"))

run_ml_pipeline(X=X[gcms_cols], y=y, model_type="random_forest",
                gcms_cols=gcms_cols,
                outfile_cv=str(output_dir / "cv_results_rf_gcms_only.pkl"),
                outfile_roc_curve=str(output_dir / "roc_rf_gcms_only.pdf"))

run_ml_pipeline(X=X[continuous_cols + binary_cols], y=y, model_type="random_forest",
                continuous_cols=continuous_cols, binary_cols=binary_cols,
                outfile_cv=str(output_dir / "cv_results_rf_clinical_only.pkl"),
                outfile_roc_curve=str(output_dir / "roc_rf_clinical_only.pdf"))

# FCNN
Xft = X.astype(np.float32)
run_ml_pipeline(X=Xft, y=y, model_type="fcnn",
                continuous_cols=continuous_cols, binary_cols=binary_cols, gcms_cols=gcms_cols,
                outfile_cv=str(output_dir / "cv_results_fcnn_gcms_and_covars.pkl"),
                outfile_roc_curve=str(output_dir / "roc_fcnn_gcms_and_covars.pdf"))

run_ml_pipeline(X=Xft[gcms_cols], y=y, model_type="fcnn",
                gcms_cols=gcms_cols,
                outfile_cv=str(output_dir / "cv_results_fcnn_gcms_only.pkl"),
                outfile_roc_curve=str(output_dir / "roc_fcnn_gcms_only.pdf"))

run_ml_pipeline(X=Xft[continuous_cols + binary_cols], y=y, model_type="fcnn",
                continuous_cols=continuous_cols, binary_cols=binary_cols,
                outfile_cv=str(output_dir / "cv_results_fcnn_clinical_only.pkl"),
                outfile_roc_curve=str(output_dir / "roc_fcnn_clinical_only.pdf"))

# ----------------------------------------
# Post-hoc: coefficient summary for lasso
# ----------------------------------------
cv_results = pd.read_pickle(output_dir / "cv_results_l1logreg_gcms_and_covars.pkl")

coefs = []
for estimator in cv_results["estimator"]:
    best_model = estimator.best_estimator_
    logreg     = best_model.named_steps["classifier"]
    coefs.append(logreg.coef_.ravel())

coefs    = np.vstack(coefs)
features = cv_results["estimator"][0].best_estimator_ \
                      .named_steps["scaler"].get_feature_names_out()

coef_summary = pd.DataFrame({
    "feature":           features,
    "mean_coef":         coefs.mean(axis=0),
    "std_coef":          coefs.std(axis=0),
    "nonzero_fraction":  (coefs != 0).mean(axis=0),
    "lower_perc12_5":    np.percentile(coefs, 12.5, axis=0),
    "upper_perc87_5":    np.percentile(coefs, 87.5, axis=0),
})
coef_summary["odds_ratio"]  = np.exp(coef_summary["mean_coef"])
coef_summary["plot_label"]  = coef_summary["feature"].str.replace("gcms__", "")
coef_summary.loc[coef_summary["nonzero_fraction"] < 0.125, "plot_label"] = ""

coef_summary.to_csv(output_dir / f"lasso_coef_summary_{experiment_name}.csv", index=False)

plt.figure(figsize=(10, 6))
sns.scatterplot(data=coef_summary, x="odds_ratio", y="nonzero_fraction",
                s=30, color="steelblue", edgecolor="black")
for _, row in coef_summary.iterrows():
    plt.text(row["odds_ratio"], row["nonzero_fraction"],
             row["plot_label"], fontsize=7, alpha=0.7)
plt.axvline(1, color="red", linestyle="--", alpha=0.6)
plt.xlabel("Odds Ratio");  plt.ylabel("Nonzero Fraction")
plt.title(f"Feature Stability — {experiment_name}")
plt.tight_layout()
plt.savefig(output_dir / f"lasso_OR_vs_nonzero_{experiment_name}.pdf",
            format="pdf", bbox_inches="tight")
plt.close()

# ----------------------------------------
# Post-hoc: Boruta for random forest
# ----------------------------------------
import numpy as np
np.int   = int
np.float = float
np.bool  = bool
from boruta import BorutaPy

cv_results = pd.read_pickle(output_dir / "cv_results_rf_gcms_and_covars.pkl")

boruta_rankings, boruta_selected, boruta_selected_weak = [], [], []

for fold_idx, estimator in enumerate(cv_results["estimator"]):
    print(f"Boruta fold {fold_idx+1}")
    best_model = estimator.best_estimator_
    rf         = best_model.named_steps["classifier"]
    train_idx  = cv_results["indices"]["train"][fold_idx]

    X_train_transformed = best_model.named_steps["scaler"].transform(X.iloc[train_idx])
    y_train_fold        = y.iloc[train_idx]

    boruta = BorutaPy(rf, n_estimators="auto", max_iter=100, random_state=42 + fold_idx)
    boruta.fit(X_train_transformed, y_train_fold.values)

    boruta_rankings.append(boruta.ranking_)
    boruta_selected.append(boruta.support_)
    boruta_selected_weak.append(boruta.support_weak_)

boruta_rankings      = np.vstack(boruta_rankings)
boruta_selected      = np.vstack(boruta_selected)
boruta_selected_weak = np.vstack(boruta_selected_weak)

features = (cv_results["estimator"][0].best_estimator_
                       .named_steps["scaler"].get_feature_names_out())

boruta_summary = pd.DataFrame({
    "feature":                 features,
    "selection_frequency":     boruta_selected.mean(axis=0),
    "weak_selection_frequency": boruta_selected_weak.mean(axis=0),
    "mean_ranking":            boruta_rankings.mean(axis=0),
    "median_ranking":          np.median(boruta_rankings, axis=0),
    "min_ranking":             boruta_rankings.min(axis=0),
    "max_ranking":             boruta_rankings.max(axis=0),
}).sort_values(["selection_frequency", "mean_ranking"], ascending=[False, True])

boruta_summary.to_csv(output_dir / f"rf_boruta_summary_{experiment_name}.csv", index=False)
print(boruta_summary.head(20))



# random forest

import pandas as pd
import numpy as np
import torch
from captum.attr import IntegratedGradients

# Load CV results
# cv_results = pd.read_pickle("/data/local/jy1008/SaMu/results/03032026/cv_results_fcnn_full.pkl")
cv_results = pd.read_pickle(output_dir / "cv_results_fcnn_gcms_and_covars.pkl")

all_importances = []

for fold_idx, estimator in enumerate(cv_results["estimator"]):
    print(f"Processing fold {fold_idx+1}/{len(cv_results['estimator'])}...")

    # Get trained pipeline
    pipeline = estimator.best_estimator_

    # Get test indices
    test_idx = cv_results["indices"]["test"][fold_idx]
    X_test_fold = X.iloc[test_idx].astype(np.float32)
    y_test_fold = y.iloc[test_idx]

    # Scale features using trained scaler
    X_input = pipeline.named_steps['scaler'].transform(X_test_fold)
    X_input = torch.tensor(X_input, dtype=torch.float32)

    # Get PyTorch model
    torch_model = pipeline.named_steps['classifier'].module_
    torch_model.eval()

    # Integrated Gradients on logits of class 1
    def logits_class1(x):
        return torch_model(x)[:, 1]  # class 1 logit

    ig = IntegratedGradients(logits_class1)
    attributions, delta = ig.attribute(
        X_input,
        target=None,
        return_convergence_delta=True
    )

    all_importances.append(attributions.detach().numpy())

# Stack fold-wise attributions
all_importances = np.vstack(all_importances)

# Feature names
features = pipeline.named_steps["scaler"].get_feature_names_out()

# Summary DataFrame
importance_summary = pd.DataFrame({
    "feature": features,
    "mean_importance": all_importances.mean(axis=0),
    "std_importance": all_importances.std(axis=0),
    "positive_fraction": (all_importances > 0).mean(axis=0),
    "lower_perc12_5": np.percentile(all_importances, 12.5, axis=0),
    "upper_perc87_5": np.percentile(all_importances, 87.5, axis=0)
})


# Save summary to CSV
importance_summary.to_csv(output_dir / f"fcnn_feature_importances_IG_{experiment_name}.csv", index=False)
print("Integrated Gradients feature importance summary saved.")





# combine results into one table

linreg_res = pd.read_csv(output_dir / f"lasso_coef_summary_{experiment_name}.csv")
linreg_res["feature"] = linreg_res["feature"].str.replace("clr_taxa__", "", regex=False)
rf_res = pd.read_csv(output_dir / f"rf_boruta_summary_{experiment_name}.csv")
rf_res["feature"] = rf_res["feature"].str.replace("clr_taxa__", "", regex=False)
fcnn_res = pd.read_csv(output_dir / f"fcnn_feature_importances_IG_{experiment_name}.csv")
fcnn_res["feature"] = fcnn_res["feature"].str.replace("clr_taxa__", "", regex=False)

merged_table = linreg_res[["feature", "mean_coef", "std_coef",
                           "odds_ratio", "nonzero_fraction",
                           "lower_perc12_5", "upper_perc87_5"]]

merged_table = merged_table.rename(columns={
    "mean_coef": "linreg_mean_coef",
    "std_coef": "linreg_std_coef",
    "odds_ratio": "linreg_odds_ratio",
    "nonzero_fraction": "linreg_nonzero_fraction",
    "lower_perc12_5": "linreg_lower_perc12_5",
    "upper_perc87_5": "linreg_upper_perc87_5"
})
merged_table = pd.merge(merged_table,
                        rf_res[["feature", "selection_frequency",
                                "weak_selection_frequency", "mean_ranking",
                                "median_ranking", "min_ranking", "max_ranking"]],
                        on="feature", how="outer")
merged_table = merged_table.rename(columns={
    "selection_frequency": "rf_selection_frequency",
    "weak_selection_frequency": "rf_weak_selection_frequency",
    "mean_ranking": "rf_mean_ranking",
    "median_ranking": "rf_median_ranking",
    "min_ranking": "rf_min_ranking",
    "max_ranking": "rf_max_ranking"
})
merged_table = pd.merge(merged_table, 
                        fcnn_res[["feature", "mean_importance",
                                  "std_importance", "positive_fraction",
                                  "lower_perc12_5", "upper_perc87_5"]],
                        on="feature", how="outer")
merged_table = merged_table.rename(columns={
    "mean_importance": "fcnn_mean_importance",
    "std_importance": "fcnn_std_importance",
    "positive_fraction": "fcnn_positive_fraction",
    "lower_perc12_5": "fcnn_lower_perc12_5",
    "upper_perc87_5": "fcnn_upper_perc87_5"
})
merged_table.to_csv(output_dir / "merged_feature_importance.csv", index=False)
