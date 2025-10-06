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
# meta_cov = meta_cov[(meta_cov['sex'] == 0)]

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

def run_lasso_logreg_pipeline(X, y, continuous_cols=None, binary_cols=None, clr_taxa_cols=None):
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
    pipeline = Pipeline([
        ('scaler', col_scaler),
        ('classifier', LogisticRegression(penalty='l1', solver='saga',
                                          max_iter=5000, random_state=42))
    ])

    # Parameter grid for inner loop to find best hyperparameters
    # NOTE: 'classifier' is the name of the pipeline step containing the logistic regression
    param_grid = {
        'classifier__C': np.logspace(-2, 2, 40),
        #'penalty': 'l1'
        #'classifier__kernel': ['linear', 'rbf']
        #'classifier__C': np.logspace(-2, 1, 30),
        #'classifier__l1_ratio': [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1.0]
    }

    # Nested cross-validation strategy
    # Outer loop: 5-fold CV
    # NOTE: in practice, I found that performance tends to vary widely based on
    # the sex balance in the k folds.
    # outer_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=111)
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=111)
    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=222)
    # outer_cv = KFold(n_splits=5, shuffle=True, random_state=101)
    # inner_cv = KFold(n_splits=5, shuffle=True, random_state=202)

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y), 1):
        y_train, y_test = y[train_idx], y[test_idx]
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

    # Plot ROC curves using RocCurveDisplay
    # https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    n_splits = outer_cv.get_n_splits()
    mean_fpr = np.linspace(0, 1, 100)
    interp_tprs = []

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    curve_kwargs_list = [dict(alpha=0.3, lw=1, color=colors[i % len(colors)]) for i in range(n_splits)]
    names = [f'ROC fold {i+1}' for i in range(n_splits)]

    mean_fpr = np.linspace(0, 1, 100)
    interp_tprs = []

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
run_lasso_logreg_pipeline(X=X, y=y, continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=taxa_cols)

# only taxonomic predictors
# X_taxon = X[taxa_cols]
run_lasso_logreg_pipeline(X=X, y=y, continuous_cols=None, binary_cols=None, clr_taxa_cols=taxa_cols)

# only clinical covariates
# X_clinical = X[clinical_cols]
run_lasso_logreg_pipeline(X=X, y=y, continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=None)



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

