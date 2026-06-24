"""loaders.py

Per-assay data loaders. Each returns a uniform contract consumed by the entry
scripts and ml_engine.run_nested_cv:

    X              : feature matrix (samples x features), covariates included
    y              : binary sarcopenia label aligned to X
    feature_groups : dict {name: [columns]} of StandardScaler-scaled groups
                     (e.g. continuous covariates, CLR taxa, metabolites)
    binary_cols    : covariate columns passed through unscaled

The engine builds its ColumnTransformer from feature_groups + binary_cols, so
adding a new assay only means writing a loader here, not touching the engine.

Shared conventions:
  * continuous covariates: age_def, nutr_score, bmi  (scaled)
  * binary covariates:     sex, smke, alco            (passthrough; sex -> 0/1)
  * young controls (age_def < 50) are dropped
  * y = 1 if sarc_status > 0 (Sarc) else 0
"""
import numpy as np
import pandas as pd

CONTINUOUS_COLS = ["age_def", "nutr_score", "bmi"]
BINARY_COLS = ["sex", "smke", "alco"]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def _clinical_covariates(meta_df, index_col):
    """Covariate frame indexed by `index_col`, with sex mapped to 0/1, young
    controls (age_def < 50) dropped, and NA covariate rows removed."""
    meta_df = meta_df[meta_df["Full.SaMu"] == 1].copy()
    meta_df["sarc_status_bin"] = (meta_df["sarc_status"] > 0).astype(int)

    cov = meta_df.set_index(index_col)[
        CONTINUOUS_COLS + BINARY_COLS + ["sarc_status_bin"]
    ].copy()
    cov = cov[~(cov["age_def"] < 50)]
    cov["sex"] = cov["sex"].map({"M": 0, "F": 1}).fillna(cov["sex"])
    return cov.dropna()


def _clr_transform(wide, pseudo_count=1e-6):
    """Centered log-ratio transform of a samples x features relative-abundance
    matrix. A pseudo-count is added and rows renormalized before the CLR."""
    x = wide + pseudo_count
    x = x.div(x.sum(axis=1), axis=0)
    geo_mean = np.exp(np.log(x).mean(axis=1))
    return np.log(x.div(geo_mean, axis=0))


def _assemble(feature_block, cov, meta_df, y_index_col, group_name):
    """Align a feature block to the covariate frame, concatenate, build y, and
    return the uniform (X, y, feature_groups, binary_cols) contract."""
    common = feature_block.index.intersection(cov.index)
    X = pd.concat([feature_block.loc[common], cov.loc[common]], axis=1)

    # y comes from the covariate frame's binarized label (already aligned to
    # `common`), avoiding a dependency on sarc_status_bin existing on meta_df.
    y = cov.loc[common, "sarc_status_bin"].astype(int)

    if "sarc_status_bin" in X.columns:
        X = X.drop(columns="sarc_status_bin")
    assert "sarc_status_bin" not in X.columns

    feature_groups = {
        "cont": CONTINUOUS_COLS,
        group_name: list(feature_block.columns),
    }
    return X, y, feature_groups, BINARY_COLS


# ---------------------------------------------------------------------------
# Metagenomics (MetaPhlAn species relative abundances)
# ---------------------------------------------------------------------------
def load_metagenomics(meta_filtered, meta_df, clr_transform=True,
                      pseudo_count=1e-6):
    """Species relative abundances pivoted wide and CLR- (or log-) transformed,
    joined to clinical covariates. Indexed by File_ID == Sample.

    With clr_transform=False the features are log(x + pseudo_count) instead,
    which is the form used for the tree-based models in the original script.
    """
    cov = _clinical_covariates(meta_df, index_col="File_ID")

    wide = meta_filtered.pivot_table(
        index="Sample", columns="Species",
        values="relative_abundance", fill_value=0,
    )
    if clr_transform:
        feats = _clr_transform(wide, pseudo_count)
    else:
        feats = np.log(wide + pseudo_count)
    assert feats.shape[0] == feats.dropna().shape[0]

    # Taxa columns are scaled as their own group ("clr_taxa") for parity with
    # the original feature-name prefixes used downstream.
    X, y, feature_groups, binary_cols = _assemble(
        feats, cov, meta_df, y_index_col="File_ID", group_name="clr_taxa")
    return X, y, feature_groups, binary_cols


# ---------------------------------------------------------------------------
# Mass-spec proteomics (DEP2 VSN-imputed wide matrix)
# ---------------------------------------------------------------------------
def load_mass_spec(data_filtered, meta_df):
    """VSN-imputed proteomics matrix (samples x proteins) joined to covariates.

    Sample IDs carry a Sarc_/NoSarc_ prefix and are integer record_ids; the
    prefix is stripped so the index matches meta_df.record_id.
    """
    cov = _clinical_covariates(meta_df, index_col="record_id")

    feats = data_filtered.copy()
    feats["Sample"] = (feats["Sample"]
                       .str.removeprefix("Sarc_")
                       .str.removeprefix("NoSarc_")
                       .astype("int64"))
    feats = feats.set_index("Sample")
    assert feats.shape[0] == feats.dropna().shape[0]

    X, y, feature_groups, binary_cols = _assemble(
        feats, cov, meta_df, y_index_col="record_id", group_name="proteins")
    return X, y, feature_groups, binary_cols


# ---------------------------------------------------------------------------
# NMR metabolites (combined with metagenomics taxa + covariates)
# ---------------------------------------------------------------------------
def load_nmr(meta_filtered, meta_df, nmr_df,
             start_col="5-Aminopentanoate", end_col="Valine",
             drop_metabolites=("Glucose", "Galactose"),
             log_transform=True, pseudo_count=1e-6):
    """CLR taxa + clinical covariates + NMR metabolites in one matrix.

    NMR metabolites are the contiguous column block [start_col, end_col],
    minus drop_metabolites, coerced numeric, NA->0, zero-sum rows dropped, then
    optionally log10(x + eps). This loader threads a third scaled group ("nmr")
    alongside the taxa and continuous-covariate groups.
    """
    cov = _clinical_covariates(meta_df, index_col="File_ID")

    wide = meta_filtered.pivot_table(
        index="Sample", columns="Species",
        values="relative_abundance", fill_value=0,
    )
    taxa = _clr_transform(wide, pseudo_count)
    assert taxa.shape[0] == taxa.dropna().shape[0]

    # Taxa + covariates first (File_ID indexed).
    common = taxa.index.intersection(cov.index)
    X = pd.concat([taxa.loc[common], cov.loc[common]], axis=1)

    # NMR block, indexed by File_ID via the meta_df join done by the caller.
    nmr = nmr_df.replace("UNK_E_1", "0").set_index("File_ID")
    cols = nmr.columns
    nmr_cols = list(cols[cols.get_loc(start_col):cols.get_loc(end_col) + 1])
    nmr_cols = [c for c in nmr_cols if c not in drop_metabolites]

    nmr_sub = nmr.loc[nmr.index.intersection(X.index), nmr_cols]
    nmr_sub = nmr_sub.apply(pd.to_numeric, errors="coerce").fillna(0)
    nmr_sub = nmr_sub[nmr_sub.sum(axis=1) > 0]
    if log_transform:
        nmr_sub = np.log10(nmr_sub + pseudo_count)

    common = X.index.intersection(nmr_sub.index)
    X = pd.concat([X.loc[common], nmr_sub.loc[common, nmr_cols]], axis=1)

    y = cov.loc[common, "sarc_status_bin"].astype(int)
    if "sarc_status_bin" in X.columns:
        X = X.drop(columns="sarc_status_bin")
    assert "sarc_status_bin" not in X.columns

    taxa_cols = [c for c in taxa.columns if c in X.columns]
    feature_groups = {
        "cont": CONTINUOUS_COLS,
        "clr_taxa": taxa_cols,
        "nmr": nmr_cols,
    }
    return X, y, feature_groups, BINARY_COLS


# ---------------------------------------------------------------------------
# LC-MS / GC-MS (DEP2 VSN-imputed wide matrix, separate metadata file)
# ---------------------------------------------------------------------------
def load_lc_gc_ms(feature_df, meta_df):
    """LC/GC-MS VSN-imputed matrix joined to its own DEP2 metadata file.

    Sample IDs are the trailing integer of the DEP2 sample label, re-prefixed
    with "SaMu" to match meta_df.label. sarc_status_bin arrives as a string
    ("Sarc"/"NoSarc"/"Unknown") in this metadata, so it is binarized here.
    Feature column names are sanitized (non-word chars -> "_").
    """
    feats = feature_df.set_index("Sample")
    feats.index = "SaMu" + feats.index.str.extract(r"(\d+)$")[0]
    feats.index.name = "Sample"

    meta = meta_df.set_index("label")
    meta = meta[meta["sarc_status_bin"].notna()]
    meta = meta[meta["sarc_status_bin"] != "Unknown"].copy()
    meta["sarc_status_bin"] = (meta["sarc_status_bin"] == "Sarc").astype(int)
    meta["sex"] = meta["sex"].map({"M": 0, "F": 1}).fillna(meta["sex"])

    common = feats.index.intersection(meta.index)
    feats, meta = feats.loc[common], meta.loc[common]

    valid = meta[CONTINUOUS_COLS + BINARY_COLS].dropna().index
    common = common.intersection(valid)
    feats, meta = feats.loc[common], meta.loc[common]

    feats.columns = feats.columns.str.replace(r"[^\w]", "_", regex=True).str.strip()
    ms_cols = list(feats.columns)

    X = pd.concat([meta[CONTINUOUS_COLS + BINARY_COLS], feats], axis=1)
    y = meta["sarc_status_bin"]
    assert "sarc_status_bin" not in X.columns
    assert X.shape[0] == y.shape[0]

    feature_groups = {"cont": CONTINUOUS_COLS, "gcms": ms_cols}
    return X, y, feature_groups, BINARY_COLS
