# SaMu sarcopenia prediction — ML models

Nested cross-validation models predicting binary sarcopenia status from each
omics assay (plus clinical covariates). The four assays share one engine; only
the data loading differs, so the per-assay entry scripts are thin.

## Layout

```
fcnn_wrapper.py        FC_NN (skorch/PyTorch) classifier + builder
ml_engine.py           Nested-CV engine: model defs, ColumnTransformer, ROC plot
feature_importance.py  Post-hoc extractors (lasso coef, Boruta RF, FC_NN IG, merge)
loaders.py             One loader per assay -> (X, y, feature_groups, binary_cols)
config.py              Per-assay paths + shared seeds/metric
run_metagenomics.py    Entry script: MetaPhlAn species
run_mass_spec.py       Entry script: DEP2 VSN proteomics
run_nmr.py             Entry script: NMR metabolites + taxa
run_lc_gc_ms.py        Entry script: LC/GC-MS features
```

## Running

Edit paths in `config.py`, then run an assay's entry script from this folder:

```bash
python run_metagenomics.py
python run_mass_spec.py
python run_nmr.py
python run_lc_gc_ms.py     # set experiment_name in config.LC_GC_MS first
```

Each writes, into that assay's `output_dir`:
`cv_results_<model>_<subset>.pkl`, `roc_curve_<model>_<subset>.pdf`, and the
feature-importance CSVs. Subsets are `full`, `<feature>_only`, `clinical_only`.

## How it fits together

A loader returns a uniform contract:

- `X` — features + clinical covariates, samples x columns
- `y` — binary label (`sarc_status > 0`)
- `feature_groups` — dict `{name: [cols]}` of StandardScaler-scaled groups
  (e.g. `cont`, `clr_taxa`, `nmr`, `proteins`, `gcms`)
- `binary_cols` — covariate columns passed through unscaled

`ml_engine.build_column_transformer` turns `feature_groups` + `binary_cols`
into the preprocessing step, so a new assay only needs a loader — the engine
never changes. `run_model_set` runs the standard full / feature-only /
clinical-only sweep per model; `run_nested_cv` is the single-model core.

Shared modelling choices (unchanged from the original scripts): outer/inner
5-fold StratifiedKFold (seeds 222 / 111), `roc_auc` scoring, the
Dawkins-et-al. lambda path for L1 logistic regression, young controls
(`age_def < 50`) dropped, `sex` mapped to 0/1.

## Per-assay quirks preserved

- **Metagenomics** — taxa pivoted wide; CLR for lasso/FC_NN, log-abundance for
  random forest (two loader calls).
- **Mass-spec** — `Sarc_`/`NoSarc_` sample-ID prefixes stripped to integer
  `record_id`.
- **NMR** — metabolite block `[5-Aminopentanoate … Valine]` minus
  `Glucose`/`Galactose`, numeric-coerced, zero-sum rows dropped, optional
  log10; threaded as a third scaled group alongside taxa.
- **LC/GC-MS** — separate DEP2 metadata file, `Unknown` diagnosis dropped,
  string `Sarc`/`NoSarc` binarized, sample IDs re-prefixed `SaMu`, feature
  names sanitized (non-word chars → `_`).

## Dependencies

numpy, pandas, scikit-learn, matplotlib, seaborn; torch + skorch (FC_NN),
captum (Integrated Gradients), boruta (RF selection). The lasso and RF paths
run without the deep-learning extras; FC_NN and IG need torch/skorch/captum,
and the Boruta summary needs boruta.

## Notes on this refactor

The four original scripts duplicated one engine around per-assay loaders. They
were consolidated into the shared `ml_engine` / `feature_importance` modules
with thin entry scripts; commented-out dead code and stale hardcoded paths were
removed. The loader-specific data handling above is preserved as-is. The lasso,
RF, ColumnTransformer, and importance-merge paths were exercised on synthetic
data; the torch/captum/boruta paths are syntax-checked but were not executed
here (those packages weren't available in the refactor environment).
