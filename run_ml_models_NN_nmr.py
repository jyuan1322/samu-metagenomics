import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, StratifiedKFold, RepeatedStratifiedKFold, GridSearchCV, cross_val_score, cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score, auc, RocCurveDisplay
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier

# for neural net
from fcnn_wrapper import *

log_transform_nmr = False

input_dir = Path("/data/local/jy1008/SaMu/results/02102026")
nmr_dir = Path("/data/local/jy1008/SaMu/metadata")

meta_filtered = pd.read_csv(input_dir / "meta_filtered.csv")
meta_df = pd.read_csv(input_dir / "meta_df_FullSaMu.csv")
nmr_df = pd.read_csv( nmr_dir / "20251027_SaMu_NMR_ureum_creat_v10.csv")

nmr_df = pd.merge(nmr_df, meta_df[['record_id', 'File_ID']], on='record_id', how='left')

# replace UNK_E_1 with 0
nmr_df = nmr_df.replace("UNK_E_1", "0")

# meta_filtered = pd.read_csv("meta_filtered.csv")
# meta_df = pd.read_csv("meta_df_FullSaMu.csv")
# nmr_df = pd.read_csv("nmr_quorum_FullSaMu.csv")

print(meta_filtered.head())
print(meta_filtered.info())
print(meta_df.head())
print(meta_df.info())

# remove samples who have Full.SaMu == 0
meta_df = meta_df[(meta_df['Full.SaMu'] == 1)]
meta_df = meta_df.copy()

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
                                         'nutr_score', 'bmi', 'stvol', 'sarc_status_bin']]
continuous_cols = ['age_def', 'nutr_score', 'bmi', 'stvol']
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
# align NMR with the other tables
nmr_df = nmr_df.set_index('File_ID')
common_samples = X.index.intersection(nmr_df.index)
# grab the metabolite columns only, and remove high-NA ones
all_cols = nmr_df.columns
# start_col = 'X5.Aminopentanoate'
start_col = '5-Aminopentanoate'
end_col = 'Valine'
start_idx = all_cols.get_loc(start_col)
end_idx = all_cols.get_loc(end_col)
nmr_cols = all_cols[start_idx:end_idx + 1]  # +1 because slicing is exclusive at the end
# to_remove = ["Glucose", "Galactose", "Urea", "Cadaverine"]
to_remove = ["Glucose", "Galactose"]
nmr_cols = [x for x in nmr_cols if x not in to_remove]
nmr_df_sub = nmr_df.loc[common_samples, nmr_cols]
# log relative abundance
# 2. Convert all columns to numeric, invalid parsing → NaN
nmr_df_sub = nmr_df_sub.apply(pd.to_numeric, errors='coerce')
# 3. Replace NaN with 0
nmr_df_sub = nmr_df_sub.fillna(0)
# remove rows whose sums are 0
nmr_df_sub = nmr_df_sub[nmr_df_sub.sum(axis=1) > 0]

# 2/11/2026 DON'T compute relative abundance
# Compute relative abundance
# row_sums = nmr_df_sub.sum(axis=1)
# nmr_relabund = nmr_df_sub.div(row_sums, axis=0)

# log-transform
if log_transform_nmr:
    eps = 1e-6
    # nmr_logrel = np.log10(nmr_relabund + eps)
    nmr_logrel = np.log10(nmr_df_sub + eps)
    common_samples = X.index.intersection(nmr_logrel.index)

    X = pd.concat([X.loc[common_samples],
                nmr_logrel.loc[common_samples, nmr_cols]], axis=1)
else:
    common_samples = X.index.intersection(nmr_df_sub.index)
    X = pd.concat([X.loc[common_samples],
                nmr_df_sub.loc[common_samples, nmr_cols]], axis=1)

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

def run_lasso_logreg_pipeline(
        X, y,
        model_type="lasso_logreg",
        continuous_cols=None,
        binary_cols=None,
        clr_taxa_cols=None,
        nmr_cols=None,
        outfile_cv="cv_results.pkl"):
    # col_scaler = ColumnTransformer(
    #     transformers=[
    #         ('cont', StandardScaler(), continuous_cols),
    #         ('bin', 'passthrough', binary_cols)
    #     ]
    # )

    # === Preprocessing ===
    transformers = []
    if continuous_cols:
        transformers.append(('cont', StandardScaler(), continuous_cols))
    if binary_cols:
        transformers.append(('bin', 'passthrough', binary_cols))
    if clr_taxa_cols:
        transformers.append(('clr_taxa', StandardScaler(), clr_taxa_cols))
    if nmr_cols:
        transformers.append(('nmr', StandardScaler(), nmr_cols))

    # print(transformers)

    # if transformers:
    #     print("applying transformers")
    #     col_scaler = ColumnTransformer(transformers=transformers)
    # else:
    #     col_scaler = 'passthrough'

    if transformers:
        col_scaler = ColumnTransformer(transformers=transformers)
        # Fit on X to compute output dimension
        X_transformed = col_scaler.fit_transform(X)
    else:
        col_scaler = 'passthrough'
        X_transformed = X

    # --- Dynamically set FC_NN input dim ---
    fcnn_input_dim = X.shape[1]
    fcnn_input_dim = X_transformed.shape[1]
    print(f"FC_NN input dimension: {fcnn_input_dim}")


    # === Model definitions ===

    # from Dawkins et al., 2024
    # https://github.com/gerberlab/mmethane/blob/main/mmethane/benchmarker.py
    lambda_min = 0.001
    path_len = 100
    l_max=100
    l_path = np.logspace(np.log10(l_max * lambda_min), np.log10(l_max), path_len)
    Cs = [1/l for l in l_path]

    # for ensemble methods
    if model_type in ['random_forest', 'adaboost']:
        grid_dict = {
            'bootstrap':[True],
            'n_estimators':[50,100],
            'max_depth':[None],
            'max_features':[None,'sqrt'],
            'min_samples_split':[2,9],
            'min_samples_leaf':[1,5],
            'learning_rate':[1e-5,5e-5,0.0001,5e-4,0.001,5e-3,0.01,0.05,0.1,0.5,1,5,10]
        }

        if model_type=='random_forest':
            rf = RandomForestClassifier()
        elif model_type=='adaboost':
            rf = AdaBoostClassifier()
        # elif model_type=='gradboost':
        #     rf = GradientBoostingClassifier()

        input_dict = rf._parameter_constraints
        rm_keys = list(set(list(grid_dict.keys())) - set(list(input_dict.keys())))
        for key in rm_keys:
            grid_dict.pop(key)

        # Prefix with classifier__
        param_grid = {f"classifier__{k}": v for k,v in grid_dict.items()}
    else:
        param_grid = None

    model_dict = {
        "lasso_logreg": {
            "estimator": LogisticRegression(
                penalty='l1', solver='saga', max_iter=5000, random_state=42
            ),
            "param_grid": {"classifier__C": Cs}
        },
        "mlp": {
            "estimator": MLPClassifier(max_iter=1000, random_state=42),
            "param_grid": {
                "classifier__hidden_layer_sizes": [(50,), (100,), (50, 50)],
                "classifier__alpha": [1e-4, 1e-3, 1e-2],
            }
        },
        "fcnn": {
            "estimator": make_ffnn_classifier(input_dim=fcnn_input_dim, h=[0, 0]),
            "param_grid": {
                "classifier__optimizer__weight_decay": [0.01],
                "classifier__module__p": [0.1],
                "classifier__lr": [0.001],
                # "classifier__lr": [1e-5],
                "classifier__max_epochs": [2000]
            }
        },
        "random_forest": {
            "estimator": RandomForestClassifier(random_state=42),
            "param_grid": param_grid
        },
        "adaboost": {
            "estimator": AdaBoostClassifier(random_state=42),
            "param_grid": param_grid
        },
    }

    if model_type not in model_dict:
        raise ValueError(f"Unknown model_type '{model_type}'. Choose from {list(model_dict.keys())}")

    model_info = model_dict[model_type]

    # === Pipeline ===
    # pipeline = Pipeline([
    #     ('scaler', col_scaler),
    #     ('classifier', LogisticRegression(penalty='l1', solver='saga',
    #                                       max_iter=5000, random_state=42))
    # ])
    pipeline = Pipeline([
        ('scaler', col_scaler),
        ('classifier', model_info["estimator"])
    ])
    param_grid = model_info["param_grid"]

    # Parameter grid for inner loop to find best hyperparameters
    # NOTE: 'classifier' is the name of the pipeline step containing the logistic regression
    # param_grid = {
    #     # 'classifier__C': np.logspace(-2, 1, 40),
    #     'classifier__C': Cs
    #     
    #     #'penalty': 'l1'
    #     #'classifier__kernel': ['linear', 'rbf']
    #     #'classifier__C': np.logspace(-2, 1, 30),
    #     #'classifier__l1_ratio': [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1.0]
    # }

    # === Cross-validation setup ===
    # Nested cross-validation strategy
    # Outer loop: 5-fold CV
    # NOTE: in practice, I found that performance improved when I set C range to [-2, 1] instead of [-3, 3]
    # outer_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=111)
    outer_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=333)

    # outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=111)
    # inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=222)
    # outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=333)
    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=444)
    # outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=101)
    # inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=202)
    # outer_cv = KFold(n_splits=5, shuffle=True, random_state=101)
    # inner_cv = KFold(n_splits=5, shuffle=True, random_state=202)

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

    # === Cross-validation ===
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
    # for neural network, disable pickling for now
    with open(outfile_cv, "wb") as f:
        pickle.dump(cv_results, f)

    # === ROC plot ===
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
# nmr cols defined above

# all predictors
run_lasso_logreg_pipeline(X=X, y=y, model_type="lasso_logreg", 
                          continuous_cols=continuous_cols,
                          binary_cols=binary_cols,
                          clr_taxa_cols=taxa_cols,
                          nmr_cols=nmr_cols,
                          outfile_cv="cv_results_l1logreg_nmr_and_taxa_full_repkfold.pkl")

# only nmr predictors and covariates
run_lasso_logreg_pipeline(X=X, y=y, model_type="lasso_logreg",
                          continuous_cols=continuous_cols,
                          binary_cols=binary_cols,
                          clr_taxa_cols=None,
                          nmr_cols=nmr_cols,
                          outfile_cv="cv_results_l1logreg_nmr_and_covars_only_repkfold.pkl")

# only nmr predictors
run_lasso_logreg_pipeline(X=X, y=y, model_type="lasso_logreg",
                          continuous_cols=None,
                          binary_cols=None,
                          clr_taxa_cols=None,
                          nmr_cols=nmr_cols,
                          outfile_cv="cv_results_l1logreg_nmr_only_repkfold.pkl")

# only clinical covariates
run_lasso_logreg_pipeline(X=X, y=y, model_type="lasso_logreg",
                          continuous_cols=continuous_cols,
                          binary_cols=binary_cols,
                          clr_taxa_cols=None,
                          nmr_cols=None,
                          outfile_cv="cv_results_l1logreg_no_nmr_clinical_only_repkfold.pkl")


# === Running random forest ===
run_lasso_logreg_pipeline(X=X, y=y, model_type="random_forest",
                          continuous_cols=continuous_cols,
                          binary_cols=binary_cols,
                          clr_taxa_cols=taxa_cols,
                          nmr_cols=nmr_cols,
                          outfile_cv="cv_results_rf_nmr_and_taxa_full_repkfold.pkl")

run_lasso_logreg_pipeline(X=X, y=y, model_type="random_forest", continuous_cols=None, binary_cols=None, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_rf_taxa_only.pkl")
run_lasso_logreg_pipeline(X=X, y=y, model_type="random_forest", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=None, outfile_cv="cv_results_rf_clinical_only.pkl")

run_lasso_logreg_pipeline(X=X, y=y, model_type="lasso_logreg", 
                          continuous_cols=continuous_cols,
                          binary_cols=binary_cols,
                          clr_taxa_cols=taxa_cols,
                          nmr_cols=nmr_cols,
                          outfile_cv="cv_results_l1logreg_nmr_and_taxa_full_repkfold.pkl")


"""
# all predictors
run_lasso_logreg_pipeline(X=X, y=y, model_type="lasso_logreg", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_l1logreg_full.pkl")

# only taxonomic predictors
# X_taxon = X[taxa_cols]
run_lasso_logreg_pipeline(X=X, y=y, model_type="lasso_logreg", continuous_cols=None, binary_cols=None, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_l1logreg_taxa_only.pkl")

# only clinical covariates
# X_clinical = X[clinical_cols]
run_lasso_logreg_pipeline(X=X, y=y, model_type="lasso_logreg", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=None, outfile_cv="cv_results_l1logreg_clinical_only.pkl")

# === Running random forest ===
run_lasso_logreg_pipeline(X=X, y=y, model_type="random_forest", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_rf_full.pkl")
run_lasso_logreg_pipeline(X=X, y=y, model_type="random_forest", continuous_cols=None, binary_cols=None, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_rf_taxa_only.pkl")
run_lasso_logreg_pipeline(X=X, y=y, model_type="random_forest", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=None, outfile_cv="cv_results_rf_clinical_only.pkl")

# === Running Adaboost ===
run_lasso_logreg_pipeline(X=X, y=y, model_type="adaboost", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_ada_full.pkl")
run_lasso_logreg_pipeline(X=X, y=y, model_type="adaboost", continuous_cols=None, binary_cols=None, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_ada_taxa_only.pkl")
run_lasso_logreg_pipeline(X=X, y=y, model_type="adaboost", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=None, outfile_cv="cv_results_ada_clinical_only.pkl")

# === Running FC_NN pytorch neural net ===
Xft = X.astype(np.float32)
run_lasso_logreg_pipeline(X=Xft, y=y, model_type="fcnn", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_fcnn_full.pkl")
run_lasso_logreg_pipeline(X=Xft, y=y, model_type="fcnn", continuous_cols=None, binary_cols=None, clr_taxa_cols=taxa_cols, outfile_cv="cv_results_fcnn_taxa_only.pkl")
run_lasso_logreg_pipeline(X=Xft, y=y, model_type="fcnn", continuous_cols=continuous_cols, binary_cols=binary_cols, clr_taxa_cols=None, outfile_cv="cv_results_fcnn_clinical_only.pkl")
"""



# extract summary stats
# cv_results = pd.read_pickle("/data/local/jy1008/SaMu/scripts/samu-metagenomics/output/10062025/cv_results_full_5fstratCV_iseed444_oseed333_Cm2_1.pkl")
# cv_results = pd.read_pickle("/data/local/jy1008/SaMu/scripts/samu-metagenomics/output/12182025/cv_results_l1logreg_nmr_and_covars_only_repkfold.pkl")
# cv_results = pd.read_pickle("/data/local/jy1008/SaMu/scripts/samu-metagenomics/output/12182025/cv_results_l1logreg_nmr_and_taxa_full_repkfold.pkl")
cv_results = pd.read_pickle("/data/local/jy1008/SaMu/results/02102026/ml_models_taxa_and_log_nmr/cv_results_l1logreg_nmr_only_repkfold.pkl")
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
# coef_summary.to_csv("lasso_L1_stratifiedkfold_i444_o333_coef_summary.csv", index=False)
# coef_summary.to_csv("lasso_L1_stratifiedkfold_nmr_and_covars_only_i444_o333_repkfold_coef_summary.csv", index=False)
# coef_summary.to_csv("lasso_L1_stratifiedkfold_nmr_and_taxa_full_i444_o333_repkfold_coef_summary.csv", index=False)
coef_summary.to_csv("lasso_L1_stratifiedkfold_nmr_only_nolog_i444_o333_repkfold_coef_summary.csv", index=False)

# run this if NMR only instead of the mask below
# coef_summary["plot_label"] = coef_summary["feature"].str.replace("nmr__", "")

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
coef_summary.loc[coef_summary["nonzero_fraction"] < 0.125, "plot_label"] = ""

"""
coef_summary["plot_label"] = coef_summary["feature"]
coef_summary.loc[coef_summary["nonzero_fraction"] < 0.125, "plot_label"] = ""
"""
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
# plt.savefig("lasso_L1_stratifiedkfold_i444_o333_OR_vs_nonzero_fraction.pdf", format="pdf", bbox_inches="tight")
# plt.savefig("lasso_L1_stratifiedkfold_nmr_and_covars_only_i444_o333_OR_vs_nonzero_fraction.pdf", format="pdf", bbox_inches="tight")
# plt.savefig("lasso_L1_stratifiedkfold_nmr_and_taxa_full_i444_o333_OR_vs_nonzero_fraction.pdf", format="pdf", bbox_inches="tight")
plt.savefig("lasso_L1_stratifiedkfold_nmr_only_nolog_i444_o333_OR_vs_nonzero_fraction.pdf", format="pdf", bbox_inches="tight")
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

