#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
spectral_analysis.py

From processed spectral XLSX → PCA/UMAP embeddings → QDA discrimination scores per dye.

Input:
    A single tidy Excel file (from spectral_preprocessing.py), with at least:
        Sample
        Amyloid
        Dye
        Excitation_nm
        Emission_nm_start
        Intensity
    Optionally:
        Intensity_norm  (normalized intensity per amyloid/dye/excitation)

Output:
    PCA and UMAP plots per dye, with QDA-based discrimination scores.

"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from collections import defaultdict

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score

import umap

# ---------------------------------------------------------------------
# 1. Column detection (simple & case-insensitive)
# ---------------------------------------------------------------------

ALIASES = {
    "amyloid": ["amyloid", "strain", "amyloid_name"],
    "dye": ["dye", "dye_name", "fluor", "fluorophore"],
    "excitation": ["excitation_nm", "excitation", "ex_nm", "ex"],
    "emission": ["emission_nm_start", "emission_nm", "emission", "lambda_em"],
    "intensity": ["intensity", "signal", "fluorescence", "value"],
    "intensity_norm": ["intensity_norm", "normalized_intensity", "norm"],
    "sample": ["sample", "sample_id", "cell", "cell_id"],
}

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
spectral_analysis.py

Input: ONE tidy-ish Excel file, containing at least these columns
(case-insensitive, names can vary slightly):

  - amyloid        e.g. "Amyloid", "strain"
  - dye            e.g. "Dye", "dye_name"
  - excitation     e.g. "Excitation_nm", "Excitation"
  - emission       e.g. "Emission_nm_start", "Emission"
  - intensity      e.g. "Intensity"

Optional:
  - normalized intensity: e.g. "Intensity_norm", "NormalizedIntensity"
  - sample ID: e.g. "Sample", "SampleID"

The script:
  - detects those columns (simple alias-based, case-insensitive)
  - runs PCA+QDA and UMAP+QDA per dye
  - shows plots colored by amyloid
"""

import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score

import umap


# ---------------------------------------------------------------------
# Column detection (simple & case-insensitive)
# ---------------------------------------------------------------------

ALIASES = {
    "amyloid": ["amyloid", "strain", "amyloid_name"],
    "dye": ["dye", "dye_name", "fluor", "fluorophore"],
    "excitation": ["excitation_nm", "excitation", "ex_nm", "ex"],
    "emission": ["emission_nm_start", "emission_nm", "emission", "lambda_em"],
    "intensity": ["intensity", "signal", "fluorescence", "value"],
    "intensity_norm": ["intensity_norm", "normalized_intensity", "norm"],
    "sample": ["sample", "sample_id", "cell", "cell_id"],
}


def find_column(df: pd.DataFrame, logical_name: str, required: bool = True):
    """
    Find a column in df that matches one of the aliases for logical_name,
    in a simple, case-insensitive way.

    Returns the actual column name, or None if not required/not found.
    """
    aliases = ALIASES[logical_name]
    lower_map = {c.lower(): c for c in df.columns}

    # exact lower-case match first
    for a in aliases:
        if a in lower_map:
            return lower_map[a]

    # fallback: substring match (very simple)
    for a in aliases:
        for lc, real in lower_map.items():
            if a in lc:
                return real

    if required:
        raise ValueError(
            f"Could not find a column for '{logical_name}'. "
            f"Looked for aliases: {aliases}. "
            f"Columns present: {list(df.columns)}"
        )
    else:
        return None


def detect_columns(df: pd.DataFrame):
    """
    Detect and return a dict of actual column names:
      cols["amyloid"], cols["dye"], cols["excitation"],
      cols["emission"], cols["intensity"],
      cols["intensity_norm"] (optional), cols["sample"] (optional)
    """
    cols = {}
    cols["amyloid"] = find_column(df, "amyloid", required=True)
    cols["dye"] = find_column(df, "dye", required=True)
    cols["excitation"] = find_column(df, "excitation", required=True)
    cols["emission"] = find_column(df, "emission", required=True)
    cols["intensity"] = find_column(df, "intensity", required=True)

    cols["intensity_norm"] = find_column(df, "intensity_norm", required=False)
    cols["sample"] = find_column(df, "sample", required=False)

    return cols


# ---------------------------------------------------------------------
# Load tidy-like Excel
# ---------------------------------------------------------------------

def load_tidy_excel(path: str):
    """Load a single Excel file and detect column names."""
    try:
        df = pd.read_excel(path, sheet_name="tidy")
    except ValueError:
        df = pd.read_excel(path)

    cols = detect_columns(df)

    # normalize amyloid/dye case for grouping
    df[cols["amyloid"]] = df[cols["amyloid"]].astype(str).str.strip().str.lower()
    df[cols["dye"]] = df[cols["dye"]].astype(str).str.strip().str.lower()

    # clean sample if present
    if cols["sample"] is not None:
        df[cols["sample"]] = df[cols["sample"]].astype(str).str.strip()

    # ------------------------------------------------------------------
    # Create new helper columns: AmyloidBase and RunNumber
    # ------------------------------------------------------------------
    amy_col = cols["amyloid"]

    # Extract base amyloid (remove trailing digits)
    df["AmyloidBase"] = (
        df[amy_col]
        .astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"\d+$", "", regex=True)
    )

    # Extract run number (digits at the end, else 1)
    df["RunNumber"] = (
        df[amy_col]
        .astype(str)
        .str.extract(r"(\d+)$")[0]
        .fillna("1")  # if no number, assume run 1
        .astype(int)
    )

    # Register in cols dictionary
    cols["amyloid_base"] = "AmyloidBase"
    cols["run_number"] = "RunNumber"
    # ------------------------------------------------------------------

    return df, cols


# ---------------------------------------------------------------------
# Build matrix per dye
# ---------------------------------------------------------------------

def build_matrix_for_dye(df: pd.DataFrame, cols: dict,
                         dye_name: str, use_normalized: bool = False):
    """
    Build feature matrix X and label vector for a single dye.

    - Filters by dye
    - Uses intensity or intensity_norm
    - Pivot: rows = samples (if sample column exists) else rows = amyloid
    - Columns = (Excitation, Emission)
    """
    amy_col = cols["amyloid"]
    amy_base_col = cols.get("amyloid_base", amy_col)   # fall back if missing
    dye_col = cols["dye"]
    ex_col = cols["excitation"]
    em_col = cols["emission"]
    int_col_raw = cols["intensity"]
    int_col_norm = cols["intensity_norm"]
    sample_col = cols["sample"]

    # choose value column
    if use_normalized and int_col_norm is not None:
        val_col = int_col_norm
    else:
        val_col = int_col_raw

    # filter by dye (case-insensitive, we already lowercased)
    dye_name = str(dye_name).lower()
    df_dye = df[df[dye_col] == dye_name].copy()
    if df_dye.empty:
        print(f"[{dye_name}] No data for this dye.")
        return None, None, False

    # index for pivot: prefer sample, else amyloid
    if sample_col is not None:
        index_key = sample_col
    else:
        index_key = amy_col

    pivot = df_dye.pivot_table(
        index=index_key,
        columns=[ex_col, em_col],
        values=val_col,
        aggfunc="mean",
    ).sort_index(axis=1)

    #imputation method
    pivot = pivot.fillna(pivot.mean())
    if pivot.shape[0] < 2:
        print(f"[{dye_name}] Not enough complete spectra after pivot/dropna.")
        return None, None, False

    X = pivot.to_numpy()

    # labels: amyloid per row
    if index_key == amy_base_col:
        labels = pivot.index.to_numpy()
    else:
        # index is sample -> map to amyloid
        idx = pivot.index
        sample_to_amy = (
            df_dye[[sample_col, amy_base_col]]
            .drop_duplicates()
            .set_index(sample_col)[amy_base_col]
        )
        labels = sample_to_amy.loc[idx].to_numpy()

    return X, labels, True


# ---------------------------------------------------------------------
# QDA helper
# ---------------------------------------------------------------------

def qda_on_embedding(X_embed, labels, n_splits: int = 10, random_state: int = 0):
    """
    Run QDA with stratified K-fold on a low-D embedding.

    Returns:
        overall_mean  : float (overall discrimination score)
        per_class_mean: dict {class: score}
    """
    X = np.asarray(X_embed)
    y = np.asarray(labels)

    classes, counts = np.unique(y, return_counts=True)
    max_splits = int(counts.min())
    if max_splits < 2:
        print("Not enough samples per class for QDA.")
        return None, None

    n_splits_eff = min(n_splits, max_splits)

    qda = QuadraticDiscriminantAnalysis()
    cv = StratifiedKFold(n_splits=n_splits_eff, shuffle=True, random_state=random_state)

    overall_acc = []
    per_class_acc = defaultdict(list)

    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        qda.fit(X_train, y_train)
        y_pred = qda.predict(X_test)

        overall_acc.append(accuracy_score(y_test, y_pred))

        for cls in classes:
            mask = (y_test == cls)
            if np.any(mask):
                per_class_acc[cls].append(
                    accuracy_score(y_test[mask], y_pred[mask])
                )

    overall_mean = float(np.mean(overall_acc))
    per_class_mean = {cls: float(np.mean(vals)) for cls, vals in per_class_acc.items()}
    return overall_mean, per_class_mean


# ---------------------------------------------------------------------
# PCA per dye
# ---------------------------------------------------------------------

def pca_for_dye(df: pd.DataFrame, cols: dict, dye_name: str,
                use_normalized: bool = False, n_components: int = 2):
    X, labels, ok = build_matrix_for_dye(df, cols, dye_name, use_normalized)
    if not ok:
        return None, None

    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X_scaled)

    overall_score, per_class = qda_on_embedding(X_pca, labels)
    if overall_score is None:
        print(f"[{dye_name}] PCA-QDA could not be computed.")
    else:
        print(f"[{dye_name}] PCA-QDA discrimination score: {overall_score:.3f}")
        for cls, s in per_class.items():
            print(f"  {cls}: {s:.3f}")

    amyloids = sorted(np.unique(labels))
    colors = plt.cm.tab10(np.linspace(0, 1, len(amyloids)))
    cmap = dict(zip(amyloids, colors))

    plt.figure(figsize=(5, 4))
    for amy in amyloids:
        mask = (labels == amy)
        plt.scatter(
            X_pca[mask, 0],
            X_pca[mask, 1],
            label=amy,
            alpha=0.5,
            c=[cmap[amy]],
        )

    var1, var2 = pca.explained_variance_ratio_[:2] * 100
    norm_tag = " (norm)" if use_normalized else ""
    title = f"PCA — dye: {dye_name}{norm_tag}"
    if overall_score is not None:
        title += f"\nQDA score: {overall_score:.3f}"
    else:
        title += f"\nQDA score: N/A"

    plt.xlabel(f"PC1 ({var1:.1f}% var)")
    plt.ylabel(f"PC2 ({var2:.1f}% var)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return overall_score, per_class


# ---------------------------------------------------------------------
# UMAP per dye
# ---------------------------------------------------------------------

def umap_for_dye(df: pd.DataFrame, cols: dict, dye_name: str,
                 use_normalized: bool = False,
                 n_neighbors: int = 15, min_dist: float = 0.1,
                 metric: str = "euclidean", random_state: int = 42):
    X, labels, ok = build_matrix_for_dye(df, cols, dye_name, use_normalized)
    if not ok:
        return None, None

    X_scaled = StandardScaler().fit_transform(X)

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
    )
    X_umap = reducer.fit_transform(X_scaled)

    overall_score, per_class = qda_on_embedding(X_umap, labels)
    if overall_score is None:
        print(f"[{dye_name}] UMAP-QDA could not be computed.")
    else:
        print(f"[{dye_name}] UMAP-QDA discrimination score: {overall_score:.3f}")
        for cls, s in per_class.items():
            print(f"  {cls}: {s:.3f}")

    amyloids = sorted(np.unique(labels))
    colors = plt.cm.tab10(np.linspace(0, 1, len(amyloids)))
    cmap = dict(zip(amyloids, colors))

    plt.figure(figsize=(5, 4))
    for amy in amyloids:
        m = labels == amy
        plt.scatter(
            X_umap[m, 0],
            X_umap[m, 1],
            label=amy,
            alpha=0.8,
            c=[cmap[amy]]
        )

    title = f"UMAP — dye: {dye_name}"
    if overall_score is not None:
        title += f"\nQDA score: {overall_score:.3f}"
    else:
        title += f"\nQDA score: N/A"

    plt.xlabel("UMAP-1")
    plt.ylabel("UMAP-2")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return overall_score, per_class

# ---------------------------------------------------------------------
# Run analysis for all dyes
# ---------------------------------------------------------------------

def run_all_dyes(tidy_path, use_normalized=False):
    df_all, cols = load_tidy_excel(tidy_path)

    dye_col = cols["dye"]   # actual column name for dye in the dataframe
    dyes = df_all[dye_col].dropna().astype(str).unique()

    pca_scores = {}
    umap_scores = {}

    for dye in sorted(dyes):
        print(f"\n=== {dye} ===")
        pca_overall, _ = pca_for_dye(df_all, cols, dye, use_normalized=use_normalized)
        umap_overall, _ = umap_for_dye(df_all, cols, dye, use_normalized=use_normalized)
        pca_scores[dye] = pca_overall
        umap_scores[dye] = umap_overall

    return df_all, pca_scores, umap_scores

def main():
    if len(sys.argv) > 1:
        tidy_path = sys.argv[1]
    else:
        tidy_path = "combined_samples_tidy.xlsx"  # change default as needed

    if not os.path.isfile(tidy_path):
        print(f"❌ Tidy-like Excel not found: {tidy_path}")
        return

    # Use normalized intensity if present?
    use_normalized = False  # set True if you want to use Intensity_norm

    run_all_dyes(tidy_path, use_normalized=use_normalized)


if __name__ == "__main__":
    main()



# # ---------------------------------------------------------------------
# # Loading processed data
# # ---------------------------------------------------------------------

# def load_processed_xlsx(pattern="*_tidy*.xlsx"):
#     """Load all processed XLSX files matching pattern and return a single df_all."""
#     files = glob(pattern)
#     if not files:
#         raise FileNotFoundError(f"No XLSX files matching pattern '{pattern}' found.")

#     dfs = []
#     for f in files:
#         df = pd.read_excel(f)
#         df["Run"] = f
#         dfs.append(df)

#     df_all = pd.concat(dfs, ignore_index=True)

#     # Clean Sample strings: "[aSyn-ThT]" -> "aSyn-ThT"
#     df_all["Sample"] = df_all["Sample"].astype(str).str.strip("[]")

#     # Split Sample into Amyloid and Dye (assuming "Amyloid-Dye" format)
#     df_all[["Amyloid", "Dye"]] = df_all["Sample"].str.split("-", n=1, expand=True)

#     return df_all

# # ---------------------------------------------------------------------
# # Matrix builder: EMBER-style feature matrix per dye
# # ---------------------------------------------------------------------

# def build_matrix_for_dye(df_all, dye_name):
#     """
#     Build feature matrix X and label vector for a single dye.

#     Returns:
#         X       : N x D numpy array
#         labels  : length-N array of amyloid names
#         emough_data      : bool, whether we had enough data
#     """
#     df_dye = df_all[df_all["Dye"] == dye_name].copy()
#     if df_dye.empty:
#         print(f"[{dye_name}] No data for this dye.")
#         return None, None, False

#     # pivot: rows = (Amyloid, Run), cols = (Excitation_nm, Emission_nm_start)
#     pivot = df_dye.pivot_table(
#         index=["Amyloid", "Run"],
#         columns=["Excitation_nm", "Emission_nm_start"],
#         values="Intensity",
#         aggfunc="mean"
#     ).sort_index(axis=1)

#     # drop rows with any missing values
#     pivot = pivot.dropna(axis=0)
#     if pivot.shape[0] < 2:
#         print(f"[{dye_name}] Not enough complete spectra after pivot/dropna.")
#         return None, None, False

#     X = pivot.to_numpy()
#     labels = pivot.index.get_level_values("Amyloid").to_numpy()
#     return X, labels, True


# """
# Quadratic Discriminant Analysis (QDA) Helper
# """

# def qda_on_embedding(X_embed, labels, n_splits=10, random_state=0):
#     """
#     Run QDA with stratified K-fold on a given 2D (or low-D) embedding.

#     Returns:
#         overall_mean  : float (overall discrimination score)
#         per_class_mean: dict {class: score}
#     """
#     X = np.asarray(X_embed)
#     y = np.asarray(labels)

#     classes, counts = np.unique(y, return_counts=True)
#     max_splits = int(counts.min())
#     if max_splits < 2:
#         print("Not enough samples per class for QDA.")
#         return None, None

#     n_splits_eff = min(n_splits, max_splits)

#     qda = QuadraticDiscriminantAnalysis()
#     cv = StratifiedKFold(n_splits=n_splits_eff, shuffle=True, random_state=random_state)

#     overall_acc = []
#     per_class_acc = defaultdict(list)

#     for train_idx, test_idx in cv.split(X, y):
#         X_train, X_test = X[train_idx], X[test_idx]
#         y_train, y_test = y[train_idx], y[test_idx]

#         qda.fit(X_train, y_train)
#         y_pred = qda.predict(X_test)

#         overall_acc.append(accuracy_score(y_test, y_pred))

#         for cls in classes:
#             mask = (y_test == cls)
#             if np.any(mask):
#                 per_class_acc[cls].append(
#                     accuracy_score(y_test[mask], y_pred[mask])
#                 )

#     overall_mean = float(np.mean(overall_acc))
#     per_class_mean = {cls: float(np.mean(vals)) for cls, vals in per_class_acc.items()}

#     return overall_mean, per_class_mean

# """
# PCA plot per dye, labeled by amyloids
# """

# def pca_for_dye(df_all, dye_name, n_components=2):
#     X, labels, enough_data = build_matrix_for_dye(df_all, dye_name)
#     if not enough_data:
#         print("not enough data")
#         return
    
#     #SCALE + PCA
#     X_scaled = StandardScaler().fit_transform(X)
#     pca = PCA(n_components=n_components)
#     X_pca = pca.fit_transform(X_scaled)

#     # QDA on the PCA coordinates
#     overall_score, per_class = qda_on_embedding(X_pca, labels)
#     if overall_score is None:
#         print(f"[{dye_name}] QDA could not be computed.")
#     else:
#         print(f"[{dye_name}] PCA-QDA discrimination score: {overall_score:.3f}")
#         for cls, s in per_class.items():
#             print(f"  {cls}: {s:.3f}")

#     #PLOT PC1 VS PC2, colored by amyloid
#     amyloids = sorted(np.unique(labels))
#     colors = plt.cm.tab10(np.linspace(0, 1, len(amyloids)))
#     cmap = dict(zip(amyloids, colors))

#     plt.figure(figsize=(5, 4))
#     for amy in amyloids:
#         mask = labels == amy
#         plt.scatter(
#             X_pca[mask, 0],
#             X_pca[mask, 1],
#             label=amy,
#             alpha=0.8,
#             c=[cmap[amy]]
#         )

#     var1, var2 = pca.explained_variance_ratio_[:2] * 100
#     title = f"PCA — dye: {dye_name}"
#     if overall_score is not None:
#         title += f"\nQDA score: {overall_score:.3f}"
#     plt.title(title)

#     plt.xlabel(f"PC1 ({var1:.1f}% var)")
#     plt.ylabel(f"PC2 ({var2:.1f}% var)")
#     plt.legend()
#     plt.tight_layout()
#     plt.show()

#     return overall_score, per_class


# """
# UMAP on full-spectrum and saw-tooth feature matrices
# - Nonlinear; preserves local neighborhood structure.
# - Often gives cleaner clusters when classes are not linearly separable.
# - Great for visualization of “strain” clusters.
# """

# def umap_for_dye(
#     df_all,
#     dye_name,
#     n_neighbors=15,
#     min_dist=0.1,
#     metric="euclidean",
#     random_state=42
# ):
#     X, labels, enough_data = build_matrix_for_dye(df_all, dye_name)
#     if not enough_data:
#         print("not enough data")
#         return

#     # often good to scale before UMAP too
#     X_scaled = StandardScaler().fit_transform(X)

#     reducer = umap.UMAP(
#         n_neighbors=n_neighbors,
#         min_dist=min_dist,
#         metric=metric,
#         random_state=random_state
#     )

#     X_umap = reducer.fit_transform(X_scaled)
#     y = np.array(labels)

#     # QDA on the UMAP coordinates
#     overall_score, per_class = qda_on_embedding(X_umap, labels)
#     if overall_score is None:
#         print(f"[{dye_name}] QDA could not be computed.")
#     else:
#         print(f"[{dye_name}] UMAP-QDA discrimination score: {overall_score:.3f}")
#         for cls, s in per_class.items():
#             print(f"  {cls}: {s:.3f}")

#     #Plot
#     amyloids = sorted(np.unique(labels))
#     colors = plt.cm.tab10(np.linspace(0, 1, len(amyloids)))
#     cmap = dict(zip(amyloids, colors))


#     plt.figure(figsize=(5, 4))
#     for amy in amyloids:
#         m = labels == amy
#         plt.scatter(
#             X_umap[m, 0],
#             X_umap[m, 1],
#             label=amy,
#             alpha=0.8,
#             c=[cmap[amy]]
#         )

#     title = f"UMAP — dye: {dye_name}"
#     if overall_score is not None:
#         title += f"\nQDA score: {overall_score:.3f}"
#     plt.xlabel("UMAP-1")
#     plt.ylabel("UMAP-2")
#     plt.title(title)
#     plt.legend()
#     plt.tight_layout()
#     plt.show()

#     return overall_score, per_class

# # ---------------------------------------------------------------------
# # Run analysis for all dyes
# # ---------------------------------------------------------------------

# def run_all_dyes(pattern="Spillover_LSM_A_processed*.xlsx"):
#     df_all = load_processed_xlsx(pattern)
#     dyes = df_all["Dye"].dropna().astype(str).unique()

#     pca_scores = {}
#     umap_scores = {}

#     for dye in sorted(dyes):
#         print(f"\n=== {dye} ===")
#         pca_overall, _ = pca_for_dye(df_all, dye)
#         umap_overall, _ = umap_for_dye(df_all, dye)
#         pca_scores[dye] = pca_overall
#         umap_scores[dye] = umap_overall

#     return df_all, pca_scores, umap_scores

# def main():
#     # Default: run PCA+UMAP+QDA on all processed files
#     run_all_dyes(pattern="Spillover_LSM_A_processed*.xlsx")

# if __name__ == "__main__":
#     main()
