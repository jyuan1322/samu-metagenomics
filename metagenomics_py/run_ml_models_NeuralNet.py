import pickle
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, StratifiedKFold, RepeatedStratifiedKFold, GridSearchCV, cross_val_score, cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score, auc, RocCurveDisplay
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay
from sklearn.neural_network import MLPClassifier

meta_filtered = pd.read_csv("meta_filtered.csv")
meta_df = pd.read_csv("meta_df.csv")

print(meta_filtered.head())
print(meta_filtered.info())
print(meta_df.head())
print(meta_df.info())

# binarize sarc_status for lasso logistic regression
meta_df["sarc_status_bin"] = (meta_df["sarc_status"] > 0).astype(int)

# Pivot to wide: rows = samples, columns = species
X_microbiome = meta_filtered.pivot_table(
    index='Sample',
    columns='Species',
    values='relative_abundance',
    fill_value=0  # species not present to 0
)
# CLR transform the relative abundances
# Replace zeros with a small pseudo-count to avoid log(0)
X_clr = X_microbiome.replace(0, 1e-6)
# Compute geometric mean per sample (row)
geometric_mean = np.exp(np.log(X_clr).mean(axis=1))
# CLR transform: log(x / gm)
X_clr = np.log(X_clr.div(geometric_mean, axis=0))

# File_ID matches Sample
# meta_cov = meta_df.set_index('File_ID')[['age_def_scaled', 'sex', 'smke', 'alco',
#                                          'nutr_score_scaled', 'bmi_scaled', 'sarc_status_bin']]
meta_cov = meta_df.set_index('File_ID')[['age_def', 'sex', 'smke', 'alco',
                                         'nutr_score', 'bmi', 'sarc_status_bin']]
continuous_cols = ['age_def', 'nutr_score', 'bmi']
binary_cols = ['sex', 'smke', 'alco']

# drop young control samples
meta_cov = meta_cov[~(meta_cov['age_def'] < 50)]

# keep only males (for testing stratification in data, and highly variable AUCs)
# meta_cov = meta_cov[(meta_cov['sex'] == 1)]

# unify samples between X and y and drop NA values
# check there are no NAs in the matrix
assert(X_clr.shape[0] == X_clr.dropna().shape[0])
# drop NA rows from covariates
meta_cov_clean = meta_cov.dropna()

# align samples
common_samples = X_clr.index.intersection(meta_cov_clean.index)
X = pd.concat([X_clr.loc[common_samples],
               meta_cov_clean.loc[common_samples]], axis=1)

y = meta_df.set_index('File_ID').loc[common_samples, 'sarc_status_bin']
y = (y > 0).astype(int)

# After NA removal, make sure sarc_status_bin is not in X
if 'sarc_status_bin' in X.columns:
    X = X.drop(columns='sarc_status_bin')

assert 'sarc_status_bin' not in X.columns

# regress X and y on the following
# - full model
# - just the rel abundances
# - just the covariates

def run_lasso_logreg_pipeline(X, y, continuous_cols=None, binary_cols=None, clr_taxa_cols=None, outfile_cv="cv_results.pkl"):
    # col_scaler = ColumnTransformer(
    #     transformers=[
    #         ('cont', StandardScaler(), continuous_cols),
    #         ('bin', 'passthrough', binary_cols)
    #     ]
    # )

    transformers = []
    if continuous_cols:
        transformers.append(('cont', StandardScaler(), continuous_cols))
    if binary_cols:
        transformers.append(('bin', 'passthrough', binary_cols))
    if clr_taxa_cols:
        transformers.append(('clr_taxa', StandardScaler(), clr_taxa_cols))

    # print(transformers)

    if transformers:
        print("applying transformers")
        col_scaler = ColumnTransformer(transformers=transformers)
    else:
        col_scaler = 'passthrough'

    # Define the model and hyperparameter grid
    """
    pipeline = Pipeline([
        ('scaler', col_scaler),
        ('classifier', LogisticRegression(penalty='l1', solver='saga',
                                          max_iter=5000, random_state=42))
    ])
    """
    pipeline = Pipeline([
    ('scaler', col_scaler),
    ('classifier', MLPClassifier(
        max_iter=5000,
        random_state=42))
    ])


    # Parameter grid for inner loop to find best hyperparameters
    # NOTE: 'classifier' is the name of the pipeline step containing the logistic regression
    """
    param_grid = {
        'classifier__C': np.logspace(-2, 1, 40),
        #'penalty': 'l1'
        #'classifier__kernel': ['linear', 'rbf']
        #'classifier__C': np.logspace(-2, 1, 30),
        #'classifier__l1_ratio': [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1.0]
    }
    """
    param_grid = {
        'classifier__hidden_layer_sizes': [(50,), (25,), (25, 25)],
        'classifier__activation': ['relu'],
        'classifier__alpha': [0.001, 0.001, 0.01, 0.1, 1],  # L2 regularization
        'classifier__learning_rate_init': [0.001, 0.01],
    }

    # Nested cross-validation strategy
    # Outer loop: 5-fold CV
    # NOTE: in practice, I found that performance improved when I set C range to [-2, 1] instead of [-3, 3]
    # outer_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=111)
    # outer_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=333)
    # outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=111)
    # inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=222)
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=333)
    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=444)
    # outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=101)
    # inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=202)
    # outer_cv = KFold(n_splits=5, shuffle=True, random_state=101)
    # inner_cv = KFold(n_splits=5, shuffle=True, random_state=202)

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y), 1):
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        print(f"Outer Fold {fold}")
        # print("  Train class distribution:", np.bincount(y_train))
        # print("  Test  class distribution:", np.bincount(y_test))
        # print(X.iloc[test_idx]["sex"].value_counts())
 
        # Train distribution
        train_table = pd.crosstab(y_train, X.iloc[train_idx]["sex"])
        print("Train distribution (y vs sex):")
        print(train_table)
        
        # Test distribution
        test_table = pd.crosstab(y_test, X.iloc[test_idx]["sex"])
        print("Test distribution (y vs sex):")
        print(test_table)
    
    print("-" * 40)

    # for fold, (train_idx, test_idx) in enumerate(inner_cv.split(X, y), 1):
    #     y_train, y_test = y[train_idx], y[test_idx]
    #     print(f"Inner Fold {fold}")
    #     print("  Train class distribution:", np.bincount(y_train))
    #     print("  Test  class distribution:", np.bincount(y_test))
    #     print(X.iloc[test_idx]["sex"].value_counts())

    scoring_metric = "roc_auc"
    # scoring_metric = "balanced_accuracy"
    # Inner loop: 5-fold CV for hyperparameter tuning
    inner_grid_search = GridSearchCV(
        estimator = pipeline,
        param_grid = param_grid,
        cv = inner_cv,
        scoring = scoring_metric,
        n_jobs = -1
    )

    # Perform nested CV
    # `cross_val_score` automatically runs the outer loop
    nested_scores = cross_val_score(
        inner_grid_search,
        X = X,
        y = y,
        cv = outer_cv,
        scoring = scoring_metric,
        n_jobs = -1
    )
    print(f"{scoring_metric}: {nested_scores}")
    # AUC scores: [0.66666667 0.74702381 0.61607143 0.55652174 0.55942029]


    # ROC curve
    cv_results = cross_validate(
        inner_grid_search,
        X=X,
        y=y,
        cv=outer_cv,
        return_estimator=True,
        return_indices=True,
        n_jobs=-1
    )
    with open(outfile_cv, "wb") as f:
        pickle.dump(cv_results, f)

    # Plot ROC curves using RocCurveDisplay
    # https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    n_splits = outer_cv.get_n_splits()
    mean_fpr = np.linspace(0, 1, 100)
    interp_tprs = []

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    curve_kwargs_list = [dict(alpha=0.3, lw=1, color=colors[i % len(colors)]) for i in range(n_splits)]
    names = [f'ROC fold {i+1}' for i in range(n_splits)]

    _, ax = plt.subplots(figsize=(6, 6))
    viz = RocCurveDisplay.from_cv_results(
        cv_results,
        X,
        y,
        ax=ax,
        name=names,
        curve_kwargs=curve_kwargs_list,
        plot_chance_level=True,
    )

    for idx in range(n_splits):
        interp_tpr = np.interp(mean_fpr, viz.fpr[idx], viz.tpr[idx])
        interp_tpr[0] = 0.0
        interp_tprs.append(interp_tpr)

    mean_tpr = np.mean(interp_tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(viz.roc_auc)

    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )

    std_tpr = np.std(interp_tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
        mean_fpr,
        tprs_lower,
        tprs_upper,
        color="grey",
        alpha=0.2,
        label=r"$\pm$ 1 std. dev.",
    )

    ax.set(
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        title=f"Mean ROC curve",
    )
    ax.legend(loc="lower right")
    plt.show()

# taxa are already CLR transformed
taxa_cols = [col for col in X.columns.tolist() if isinstance(col, str) and col.startswith("k__Bacteria")]
clinical_cols = [col for col in X.columns.tolist() if isinstance(col, str) and not col.startswith("k__Bacteria")]

# all predictors
run_lasso_logreg_pipeline(X=X, y=y, continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_full.pkl")

# only taxonomic predictors
# X_taxon = X[taxa_cols]
run_lasso_logreg_pipeline(X=X, y=y, continuous_cols=None, binary_cols=None, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_taxa_only.pkl")

# only clinical covariates
# X_clinical = X[clinical_cols]
run_lasso_logreg_pipeline(X=X, y=y, continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=None, outfile_cv="cv_results_clinical_only.pkl")



# extract summary stats
cv_results = pd.read_pickle("/data/local/jy1008/SaMu/scripts/samu-metagenomics/output/10062025/cv_results_full_5fstratCV_iseed444_oseed333_Cm2_1.pkl")
for i, estimator in enumerate(cv_results["estimator"], 1):
    best_model = estimator.best_estimator_  # the pipeline with tuned C
    logreg = best_model.named_steps["classifier"]
    coef = logreg.coef_.ravel()  # 1D array
    features = best_model.named_steps["scaler"].get_feature_names_out()
    coef_df = pd.DataFrame({"feature": features, "coef": coef})
    print(f"\n=== Outer Fold {i} ===")
    print(coef_df[coef_df["coef"] != 0].sort_values("coef", ascending=False))

coefs = []
for estimator in cv_results["estimator"]:
    best_model = estimator.best_estimator_
    logreg = best_model.named_steps["classifier"]
    coefs.append(logreg.coef_.ravel())

coefs = np.vstack(coefs)

# Just take feature names from the first model
features = cv_results["estimator"][0].best_estimator_.named_steps["scaler"].get_feature_names_out()

coef_summary = pd.DataFrame({
    "feature": features,
    "mean_coef": coefs.mean(axis=0),
    "std_coef": coefs.std(axis=0),
    "nonzero_fraction": (coefs != 0).mean(axis=0)
})
coef_summary["odds_ratio"] = np.exp(coef_summary["mean_coef"])
coef_summary.sort_values("nonzero_fraction", ascending=False).head(20)
coef_summary.to_csv("lasso_L1_stratifiedkfold_i444_o333_coef_summary.csv", index=False)

# Extract genus and species into label for taxonomic rows
mask = coef_summary["feature"].str.contains(r"g__.*\|s__", na=False)
# Initialize column as copy (so unaffected features stay intact)
coef_summary["genus_species"] = coef_summary["feature"]
coef_summary.loc[mask, "genus_species"] = (
    coef_summary.loc[mask, "feature"]
    .str.extract(r"g__([^|]+)\|s__([^|]+)")
    .fillna("")
    .agg(" ".join, axis=1)
    .str.strip()
)
coef_summary["plot_label"] = coef_summary["genus_species"]
coef_summary.loc[coef_summary["nonzero_fraction"] < 0.2, "plot_label"] = ""
# top_features = coef_summary.loc[coef_summary.nonzero_fraction > 0.3]
# top_features = top_features.sort_values("mean_coef")

# plt.figure(figsize=(6, len(top_features)/2))
# plt.barh(top_features["feature"], top_features["mean_coef"])
# plt.xlabel("Mean log-odds coefficient")
# plt.title("Stable features across folds")
# plt.tight_layout()
# plt.show()

plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=coef_summary,
    x="odds_ratio",
    y="nonzero_fraction",
    s=30,
    color="steelblue",
    edgecolor="black"
)

# Label points
for _, row in coef_summary.iterrows():
    plt.text(
        row["odds_ratio"],
        row["nonzero_fraction"],
        row["plot_label"],
        fontsize=7,
        alpha=0.7
    )

plt.axvline(1, color="red", linestyle="--", alpha=0.6)
plt.xlabel("Odds Ratio")
plt.ylabel("Nonzero Fraction")
plt.title("Feature Nonzero Fraction vs Odds Ratio")
plt.tight_layout()
plt.savefig("lasso_L1_stratifiedkfold_i444_o333_OR_vs_nonzero_fraction.pdf", format="pdf", bbox_inches="tight")
plt.show()






# Read the DESeq2 results file
# deseq2_df = pd.read_csv("output/deseq2_results.csv", index_col=0)
deseq2_df = pd.read_csv("/data/local/jy1008/SaMu/scripts/samu-metagenomics/output/10062025/deseq2_results_sarc_bin.csv", index_col=0)

# Read L1 logistic regression coefficients
# (from your coef_summary DataFrame)
logreg = coef_summary.copy()
logreg = logreg.set_index("feature")
logreg.index = logreg.index.str.replace('^clr_taxa__', '', regex=True)

merged = logreg.merge(
    deseq2_df,
    left_index=True,
    right_index=True,
    how="inner",
    suffixes=("_logreg", "_deseq")
)
print(merged.shape)

# compare direction and magnitude of effects
# Filter merged DataFrame for significant DESeq2 results
sig_merged = merged[merged["padj"] < 0.05]
sig_merged = sig_merged[sig_merged["nonzero_fraction"] > 0]

# Data
x = np.sign(sig_merged["mean_coef"].values) * sig_merged["nonzero_fraction"].values
y = sig_merged["log2FoldChange"].values
labels = sig_merged["genus_species"].values

# Scatter plot
plt.figure(figsize=(8,6))
plt.scatter(x, y, alpha=0.7)

# Trendline through origin
m = np.sum(x * y) / np.sum(x**2)
y_pred = m * x
plt.plot(x, y_pred, color='red', linewidth=2, label=f"y = {m:.2f}x")

# R² (through origin)
ss_res = np.sum((y - y_pred)**2)
ss_tot = np.sum(y**2)
r2 = 1 - ss_res/ss_tot
plt.text(0.05, 0.95, f"y = {m:.2f}x\nR² = {r2:.2f}",
         transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7))

# Add labels to points
for i, point_label in enumerate(sig_merged["genus_species"]):
  plt.text(x[i], y[i], point_label, fontsize=8, alpha=0.7)

# Axes lines
plt.axhline(0, color='gray', linestyle='--')
plt.axvline(0, color='gray', linestyle='--')
plt.xlabel("L1 Logistic sign(mean coef) * nonzero fraction")
plt.ylabel("DESeq2 log2 fold change")
plt.title("Comparison of Logistic Regression and DESeq2 effects")

plt.legend()
plt.tight_layout()
plt.savefig("lasso_L1_stratifiedkfold_i444_o333_vs_deseq2_effects.pdf", format="pdf", bbox_inches="tight")
plt.show()

sig_merged.to_csv("lasso_L1_stratifiedkfold_i444_vs_deseq2_effects.csv", index=False)

"""clf = LogisticRegressionCV(
    Cs=10,            # number of regularization strengths to try
    cv=5,             # 5-fold cross-validation
    penalty='l1',     # Lasso
    solver='saga',    # saga supports l1 penalty
    scoring='roc_auc',# optimize AUC
    max_iter=5000,
    refit=True,
    random_state=42
)
clf.fit(X_clean, y_clean)


# Predict probabilities
y_prob = clf.predict_proba(X_clean)[:, 1]

# Compute AUC
roc_auc = roc_auc_score(y_clean, y_prob)
print("AUC:", roc_auc) # AUC: 0.9112523540489642

# ROC curve
fpr, tpr, thresholds = roc_curve(y_clean, y_prob)

plt.figure(figsize=(6,6))
plt.plot(fpr, tpr, color='blue', label=f'AUC = {roc_auc:.2f}')
plt.plot([0,1], [0,1], color='grey', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc='lower right')
# plt.show()
plt.savefig("samu_lasso_logisticregression_roc_curve.pdf", format="pdf", bbox_inches="tight")
plt.close()

# run with only genomic info

# LogisticRegressionCV doesn't handle NA values
# Combine features + target
X_clr['sarc_status_bin'] = y

# Drop rows with any NaN
X_clr_clean = X_clr.dropna()
y_clr_clean = X_clr_clean.pop('sarc_status_bin')

clf = LogisticRegressionCV(
    Cs=10,            # number of regularization strengths to try
    cv=5,             # 5-fold cross-validation
    penalty='l1',     # Lasso
    solver='saga',    # saga supports l1 penalty
    scoring='roc_auc',# optimize AUC
    max_iter=5000,
    refit=True,
    random_state=42
)
clf.fit(X_clr_clean, y_clr_clean)


# Predict probabilities   
y_clr_prob = clf.predict_proba(X_clr_clean)[:, 1]

# Compute AUC
roc_auc = roc_auc_score(y_clr_clean, y_clr_prob)
print("AUC:", roc_auc) # AUC: 0.9187344913151365

# ROC curve
fpr, tpr, thresholds = roc_curve(y_clr_clean, y_clr_prob)

plt.figure(figsize=(6,6))
plt.plot(fpr, tpr, color='blue', label=f'AUC = {roc_auc:.2f}')
plt.plot([0,1], [0,1], color='grey', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc='lower right')
# plt.show()
plt.savefig("samu_lasso_logisticregression_roc_curve_genome_only.pdf", format="pdf", bbox_inches="tight")
plt.close()"""

