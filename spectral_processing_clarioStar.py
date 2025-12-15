#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
spectral_processing_clarioStar.py

Process BMG ClarioStar plate-reader *ex.xlsx emission files into:
  1) EMBER-style concatenated “sawtooth” spectra per sample
  2) Sawtooth plots (PNG) for each sample
  3) A tidy Excel file with long-format spectra

This is the plate-reader analogue of spectral_preprocessing.py
for your micrograph EMBER data.

Expected input
--------------
- One or more Excel files like: 350ex.xlsx, 405ex.xlsx, 470ex.xlsx, ...
- Each file has:
    Column 0: "Well" / "Content" (contains "Raw Data  (Em Spectrum)")
    Column 1: "Wavelength [nm]"
    Columns 2+: well IDs like "A01", "B01", ... "F04"

We:
  * auto-detect well columns (A01, B02, ...),
  * parse them into (dye letters, numeric position),
  * group by numeric position (01, 02, ...) as biological samples,
  * hard-code mapping of position → sample name below.

Outputs
-------
1. sawtooth_plots/   (folder of per-sample PNG plots)
2. clario_tidy.xlsx  (tidy long-format data)
3. clario_ember_matrix.csv (optional: wide matrix, one row per sample)

project_folder/
│
├── 350ex.xlsx
├── 405ex.xlsx
├── 470ex.xlsx
│
├── spectral_processing_clarioStar.py
│
└── clario_outputs/
     ├── clario_tidy.xlsx
     ├── clario_ember_matrix.csv
     └── sawtooth_plots/
          ├── sawtooth_Beads.png
          ├── sawtooth_Tau_fibrils.png
          └── ...

          
Edit ONLY the "USER CONFIGURATION" section to adapt to new experiments.
"""

import os
import sys
import re
import shutil
from typing import Dict, List, Tuple
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# USER CONFIGURATION (EDIT THIS BLOCK FOR NEW EXPERIMENTS)
# -------------------------------------------------------------------

# Where to look for *ex.xlsx files; "." means current directory
DATA_DIR = "."

# Central export directory for all outputs
EXPORT_DIR = f"clario_outputs_{datetime.now():%Y%m%d_%H%M%S}"

# Folder names and filenames relative to EXPORT_DIR
#SAWTOOTH_DIR = os.path.join(EXPORT_DIR, "sawtooth_plots")
SAWTOOTH_GRID_DIR = os.path.join(EXPORT_DIR, "sawtooth_grids")
SAWTOOTH_GRID_DIR_NORM = os.path.join(EXPORT_DIR, "sawtooth_grids_norm")
SAWTOOTH_GRID_DIR_BG = os.path.join(EXPORT_DIR, "sawtooth_grids_bgsub")
SAWTOOTH_GRID_DIR_BG_NORM = os.path.join(EXPORT_DIR, "sawtooth_grids_bgsub_norm")
TIDY_XLSX = os.path.join(EXPORT_DIR, "clario_tidy.xlsx")
TIDY_NORM_XLSX = os.path.join(EXPORT_DIR, "clario_tidy_norm.xlsx")
EMBER_MATRIX_CSV = os.path.join(EXPORT_DIR, "clario_ember_matrix.csv")

# Hard-coded sample names in the order of plate "positions" (01, 02, 03, ...)
# If you have 4 different amyloids in columns 01–04, for example:
SAMPLE_NAMES_IN_POSITION_ORDER = [
    "Beads",
    "Monomeric Tau",
    "Tau fibrils",
    "a-Syn fibrils",
] 

# Row identifiers used to detect spectrum rows in the first column
RAW_SUBSTR_1 = "raw data"
RAW_SUBSTR_2 = "em spectrum"

# Regex for well IDs like A01, B12, AA03 etc.
WELL_PATTERN = re.compile(r"^[A-Za-z]+(\d+)$")

# Regex to get excitation wavelength from "350ex.xlsx"
EXCITATION_PATTERN = re.compile(r"(\d+)ex\.xlsx$", re.IGNORECASE)

# -------------------------------------------------------------------
# LOW-LEVEL HELPERS
# -------------------------------------------------------------------

def extract_excitation_from_filename(fname: str):
    m = EXCITATION_PATTERN.search(fname)
    return int(m.group(1)) if m else None


def find_excel_files(data_dir: str) -> List[str]:
    """Find all *ex.xlsx files under DATA_DIR."""
    paths = []
    for root, dirs, files in os.walk(data_dir):
        for fn in files:
            if fn.lower().endswith(".xlsx") and "ex" in fn.lower():
                paths.append(os.path.join(root, fn))
    # unique and sorted by excitation wavelength
    paths = sorted(
        set(paths),
        key=lambda p: extract_excitation_from_filename(os.path.basename(p)) or 99999,
    )
    return paths


def detect_structure(df: pd.DataFrame):
    """
    Detect structure of a ClarioStar *ex.xlsx file.

    Assumes:
      - Column 0: 'Well'/'Content' (row labels, e.g. 'Raw Data  (Em Spectrum)')
      - Column 1: emission wavelength (nm)
      - Columns 2+: well IDs like 'A01', 'B01', etc.

    Returns:
      content_col    : column label for content / well
      wavelength_col : column label for emission wavelength
      well_cols      : list of well columns
      mapping        : {col_name: (dye_letters, position_string)}
      dyes           : sorted list of dye IDs
      positions      : sorted list of position strings (e.g. '01','02')
    """
    cols = list(df.columns)

    if len(cols) < 3:
        raise ValueError("Expected at least 3 columns (Well, Wavelength, wells...).")

    # Force these assignments:
    content_col = cols[0]      # "Well" or "Content"
    wavelength_col = cols[1]   # "Wavelength [nm]" or similar

    well_cols = []
    mapping = {}
    dyes_set = set()
    positions_set = set()

    # All remaining columns are potential wells like A01, B02, ...
    for c in cols[2:]:
        m = WELL_PATTERN.match(str(c))
        if not m:
            continue
        pos = m.group(1)             # numeric part
        dye = str(c)[:-len(pos)]     # letter part
        well_cols.append(c)
        mapping[c] = (dye, pos)
        dyes_set.add(dye)
        positions_set.add(pos)

    if not well_cols:
        raise ValueError("No well columns matching pattern like 'A01' were found.")

    dyes = sorted(dyes_set)
    positions = sorted(positions_set)
    return content_col, wavelength_col, well_cols, mapping, dyes, positions


def read_plate_file(path: str,
                    expected_dyes=None,
                    expected_positions=None) -> Dict:
    """
    Read one *ex.xlsx plate file and return dictionary describing it.

    Result keys:
      'path'       : original path
      'excitation' : int or None
      'waves'      : np.array of emission wavelengths
      'blocks'     : {position: [array for dye1, array for dye2, ...]}
      'dyes'       : list of dye names, in order
      'positions'  : list of position ids, in order
    """
    df = pd.read_excel(path)
    content_col, wav_col, _, mapping, dyes, positions = detect_structure(df)

    if expected_dyes is not None:
        dyes = expected_dyes
    if expected_positions is not None:
        positions = expected_positions

    content_str = df[content_col].astype(str).str.lower()
    mask = content_str.str.contains(RAW_SUBSTR_1) & content_str.str.contains(RAW_SUBSTR_2)
    if not mask.any():
        # fall back to anything containing 'raw' if more permissive needed
        mask = content_str.str.contains("raw")

    spec = df[mask]
    if spec.empty:
        raise ValueError(f"No spectrum rows detected in {path}")

    waves = spec[wav_col].astype(float).to_numpy()
    n_wave = len(waves)

    # map (dye, position) → column
    dp_to_col = {}
    for col, (dye, pos) in mapping.items():
        dp_to_col[(dye, pos)] = col

    blocks = {pos: [] for pos in positions}
    for pos in positions:
        for dye in dyes:
            col = dp_to_col.get((dye, pos), None)
            if col is not None and col in spec.columns:
                vals = spec[col].astype(float).to_numpy()
            else:
                vals = np.full(n_wave, np.nan)
            blocks[pos].append(vals)

    exc = extract_excitation_from_filename(os.path.basename(path))
    print(f"   First few emission wavelengths (nm): {waves[:5]}")
    print(f"   Example intensity A01 (first 5): {blocks[positions[0]][0][:5]}")

    return {
        "path": path,
        "excitation": exc,
        "waves": waves,
        "blocks": blocks,
        "dyes": dyes,
        "positions": positions,
    }

# -------------------------------------------------------------------
# BUILD EMBER MATRIX + SAWTOOTH + TIDY
# -------------------------------------------------------------------

def assign_sample_names(positions: List[str]) -> Dict[str, str]:
    """
    Map numeric positions ('01','02',...) to human-readable sample names
    using SAMPLE_NAMES_IN_POSITION_ORDER. Extra positions get default names.
    """
    mapping = {}
    for i, pos in enumerate(positions):
        if i < len(SAMPLE_NAMES_IN_POSITION_ORDER):
            mapping[pos] = SAMPLE_NAMES_IN_POSITION_ORDER[i]
        else:
            mapping[pos] = f"Sample {pos}"
    return mapping


def build_ember_matrix(file_data_list: List[Dict],
                       positions: List[str],
                       dyes: List[str],
                       pos_to_name: Dict[str, str]):
    """
    From all plate files, build:
      X      : (n_samples, n_features) EMBER matrix (min-max normalized)
      sample_labels : list of sample names (one per row of X)
      waves_per_ex  : list of wavelength arrays (one per excitation file)
      excitations   : list of excitation wavelengths (same order as waves_per_ex)
    """
    waves_per_ex = [fd["waves"] for fd in file_data_list]
    excitations = [fd["excitation"] for fd in file_data_list]

    sample_labels = []
    profiles = []

    for pos in positions:
        name = pos_to_name[pos]
        blocks_for_sample = []
        for fd in file_data_list:
            blocks_for_sample.extend(fd["blocks"][pos])  # dye-major
        vec = np.concatenate(blocks_for_sample)
        sample_labels.append(name)
        profiles.append(vec)

    X = np.vstack(profiles)

    # per-sample min–max normalization
    X_min = X.min(axis=1, keepdims=True)
    X_max = X.max(axis=1, keepdims=True)
    X_norm = (X - X_min) / (X_max - X_min + 1e-9)

    return X_norm, sample_labels, waves_per_ex, excitations


def ensure_dir(path: str):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def plot_sawtooth_grid_for_dye(tidy: pd.DataFrame,
                               dye_name: str,
                               use_norm: bool,
                               out_dir: str,
                               sample_order: List[str] = None):
    """
    Make a grid of spectra for one dye:

        columns = excitation wavelengths
        rows    = samples (amyloids)

    Each panel is emission_nm vs Intensity or Intensity_norm.
    """

    df_d = tidy[tidy["Dye"] == dye_name].copy()
    if df_d.empty:
        return

    # Unique excitations and samples
    ex_values = sorted(df_d["Excitation_nm"].dropna().unique())
    if sample_order is None:
        sample_values = sorted(df_d["Sample"].dropna().unique())
    else:
        # keep only samples actually present
        sample_values = [s for s in sample_order if s in df_d["Sample"].unique()]

    n_rows = len(sample_values)
    n_cols = len(ex_values)
    if n_rows == 0 or n_cols == 0:
        return

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4 * n_cols, 2.5 * n_rows),
        sharex=True,
        sharey=use_norm,   # for normalized plots, share y-axis
    )

    # In case of single row/col, axes might not be 2D
    if n_rows == 1:
        axes = np.expand_dims(axes, axis=0)
    if n_cols == 1:
        axes = np.expand_dims(axes, axis=1)

    value_col = "Intensity_norm" if use_norm else "Intensity"
    title_suffix = " (normalized)" if use_norm else ""

    for r, sample in enumerate(sample_values):
        for c, exc in enumerate(ex_values):
            ax = axes[r, c]
            sub = df_d[(df_d["Sample"] == sample) &
                       (df_d["Excitation_nm"] == exc)].copy()
            if sub.empty:
                ax.set_visible(False)
                continue

            sub = sub.sort_values("Emission_nm")
            x = sub["Emission_nm"].to_numpy()
            y = sub[value_col].to_numpy()

            ax.plot(x, y)

            if r == 0:
                ax.set_title(f"{exc} nm")
            if c == 0:
                ylabel = sample
                if use_norm:
                    ylabel += "\n(norm. max = 1)"
                ax.set_ylabel(ylabel)

            if r == n_rows - 1:
                ax.set_xlabel("Emission (nm)")

    fig.suptitle(f"Current Sawtooths – {dye_name}{title_suffix}", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    ensure_dir(out_dir)
    safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", dye_name)
    norm_tag = "_norm" if use_norm else ""
    out_path = os.path.join(out_dir, f"sawtooths_{safe_name}{norm_tag}.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)



def make_tidy_dataframe(file_data_list: List[Dict],
                        positions: List[str],
                        plate_dyes: List[str],
                        dye_labels: List[str],
                        pos_to_name: Dict[str, str]) -> pd.DataFrame:
    """
    Build a long-format DataFrame with columns:
      Sample, Position, Dye, Excitation_nm, Emission_nm, Intensity, SourceFile
    """
    records = []
    for fd in file_data_list:
        exc = fd["excitation"]
        waves = fd["waves"]
        n_wave = len(waves)
        for pos in positions:
            sample_name = pos_to_name[pos]
            blocks = fd["blocks"][pos]  # length == len(plate_dyes)
            for dye_idx, _ in enumerate(plate_dyes):
                intensities = blocks[dye_idx]
                dye_name = dye_labels[dye_idx]
                for i in range(n_wave):
                    records.append({
                        "Sample": sample_name,
                        "Position": pos,
                        "Dye": dye_name,
                        "Excitation_nm": exc,
                        "Emission_nm": waves[i],
                        "Intensity": intensities[i],
                        "SourceFile": os.path.basename(fd["path"]),
                    })
    tidy = pd.DataFrame.from_records(records)
    return tidy

def add_background_subtraction(
    tidy: pd.DataFrame,
    background_sample: str = "Beads"
) -> pd.DataFrame:
    """
    Add background-subtracted intensity columns using a reference sample.

    For each (Dye, Excitation_nm, Emission_nm) triple, we take the intensity
    from `background_sample` (e.g. "Beads") as background and subtract it from
    all samples at that same triple.

      Intensity_bgsub      = Intensity - Intensity(background)
      Intensity_bgsub_norm = per (Sample, Dye, Excitation_nm) max-normalized
                             version of Intensity_bgsub.
    """
    tidy = tidy.copy()

    key_cols = ["Dye", "Excitation_nm", "Emission_nm"]

    # Background series: intensity for the background sample at each (Dye, Exc, Em)
    bg = (
        tidy[tidy["Sample"] == background_sample]
        .set_index(key_cols)["Intensity"]
        .rename("Background_Intensity")
    )

    # Attach background intensity to every row (left-join on key_cols)
    tidy = tidy.merge(
        bg.reset_index(),
        on=key_cols,
        how="left",
    )

    # If no background found for some triple, assume 0
    tidy["Background_Intensity"] = tidy["Background_Intensity"].fillna(0.0)

    # Background-subtracted raw intensity
    tidy["Intensity_bgsub"] = tidy["Intensity"] - tidy["Background_Intensity"]

    # Now normalized version of the background-subtracted intensity
    def _norm(x):
        m = x.max()
        return x / m if m > 0 else 0.0

    tidy["Intensity_bgsub_norm"] = (
        tidy.groupby(["Sample", "Dye", "Excitation_nm"])["Intensity_bgsub"]
            .transform(_norm)
    )

    return tidy


def add_normalized_intensity(tidy: pd.DataFrame) -> pd.DataFrame:
    """
    Add Intensity_norm column:
      for each (Sample, Dye, Excitation_nm) group,
      divide intensities by the max in that group.
    Mirrors spectral_preprocessing.py (which used Amyloid).
    """
    tidy = tidy.copy()

    def _norm(x):
        m = x.max()
        return x / m if m > 0 else 0.0

    tidy["Intensity_norm"] = (
        tidy.groupby(["Sample", "Dye", "Excitation_nm"])["Intensity"]
            .transform(_norm)
    )
    return tidy


# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

def main():

    print("=" * 72)
    print("  ClarioStar spectral processing → EMBER + sawtooth + tidy XLSX")
    print("=" * 72)

    # Create export directories
    if os.path.exists(EXPORT_DIR):
        shutil.rmtree(EXPORT_DIR)

    ensure_dir(EXPORT_DIR)
    ensure_dir(SAWTOOTH_GRID_DIR)
    ensure_dir(SAWTOOTH_GRID_DIR_NORM)
    ensure_dir(SAWTOOTH_GRID_DIR_BG)
    ensure_dir(SAWTOOTH_GRID_DIR_BG_NORM)


    excel_paths = find_excel_files(DATA_DIR)
    if not excel_paths:
        print("\n❌ No *ex.xlsx files found under", os.path.abspath(DATA_DIR))
        return

    print("\n📊 Found Excel files:")
    for p in excel_paths:
        print(f"   {os.path.basename(p):20s}  →  "
              f"{extract_excitation_from_filename(os.path.basename(p))} nm")

    # Read first file to infer structure (positions + instrument dye codes)
    first_df = pd.read_excel(excel_paths[0])
    _, _, _, _, plate_dyes, positions = detect_structure(first_df)

    # Human-readable dye names in the same order as the plate rows (A, B, C, ...)
    dye_labels = ["NIAD4", "LDS698", "BF188", "DQTCI", "MCAAD3", "ThT"]

    if len(dye_labels) != len(plate_dyes):
        raise ValueError(
            f"Number of dye labels ({len(dye_labels)}) does not match "
            f"number of plate dye codes ({len(plate_dyes)}: {plate_dyes})"
        )

    print("\nUsing dye mapping:")
    print("  Plate dyes:   ", ", ".join(plate_dyes))
    print("  Dye labels:   ", ", ".join(dye_labels))
    print("  Positions:    ", ", ".join(positions))


    # Hard-coded sample names (no terminal prompts)
    pos_to_name = assign_sample_names(positions)
    print("\nSample name mapping (position → name):")
    for pos in positions:
        print(f"  {pos:>2s} → {pos_to_name[pos]}")

    # Read all plate files with consistent dye / position order
    file_data_list = []
    for p in excel_paths:
        print(f"\n📄 Reading {p}")
        try:
            fd = read_plate_file(p, expected_dyes=plate_dyes, expected_positions=positions)
            file_data_list.append(fd)
            print(f"   ✓ Excitation {fd['excitation']} nm, "
                  f"{len(fd['waves'])} emission points")
        except Exception as e:
            print(f"   ✗ Error: {e}")

    if not file_data_list:
        print("\n❌ No valid plate files were read.")
        return

    # Build EMBER matrix
    X_ember, sample_labels, waves_per_ex, excitations = build_ember_matrix(
        file_data_list, positions, plate_dyes, pos_to_name
    )

    # Save EMBER matrix (optional)
    if EMBER_MATRIX_CSV:
        ember_df = pd.DataFrame(
            X_ember,
            index=sample_labels,
            columns=[f"feat_{i+1}" for i in range(X_ember.shape[1])]
        )
        ember_df.to_csv(EMBER_MATRIX_CSV)
        print(f"\n💾 EMBER matrix saved to {EMBER_MATRIX_CSV}")

    # ---------------------------------------------------------
    # Build tidy dataframe (long format) and normalized version
    # ---------------------------------------------------------
    print("\n📑 Building tidy dataframe…")
    tidy = make_tidy_dataframe(
        file_data_list, positions, plate_dyes, dye_labels, pos_to_name
    )

    # Original normalization (per Sample × Dye × Exc)
    tidy = add_normalized_intensity(tidy)

    # Beads-subtracted intensities + their own normalization
    tidy = add_background_subtraction(tidy, background_sample="Beads")

    # Save one tidy file with all columns
    tidy.to_excel(TIDY_XLSX, index=False)
    print(f"   ✓ Tidy spectra (incl. Intensity_norm, Intensity_bgsub, Intensity_bgsub_norm) saved to {TIDY_XLSX}")

    # A copy of tidy where 'Intensity'/'Intensity_norm' are the BEADS-SUBTRACTED values
    tidy_bg = tidy.copy()
    tidy_bg["Intensity"] = tidy_bg["Intensity_bgsub"]
    tidy_bg["Intensity_norm"] = tidy_bg["Intensity_bgsub_norm"]


    # ---------------------------------------------------------
    # Generate grid-style "current sawtooth" plots by dye
    # ---------------------------------------------------------
    print("\n📈 Generating grid sawtooth plots by dye…")

    # Use the same sample ordering as in SAMPLE_NAMES_IN_POSITION_ORDER where possible
    sample_order = SAMPLE_NAMES_IN_POSITION_ORDER

    for dye_name in dye_labels:
        print(f"   Dye: {dye_name}")

        # 1) Raw curves (no background subtraction)
        plot_sawtooth_grid_for_dye(
            tidy=tidy,
            dye_name=dye_name,
            use_norm=False,
            out_dir=SAWTOOTH_GRID_DIR,
            sample_order=sample_order,
        )

        # 2) Normalized curves (no background subtraction)
        plot_sawtooth_grid_for_dye(
            tidy=tidy,
            dye_name=dye_name,
            use_norm=True,
            out_dir=SAWTOOTH_GRID_DIR_NORM,
            sample_order=sample_order,
        )

        # 3) Beads-subtracted raw curves
        plot_sawtooth_grid_for_dye(
            tidy=tidy_bg,
            dye_name=dye_name,
            use_norm=False,
            out_dir=SAWTOOTH_GRID_DIR_BG,
            sample_order=sample_order,
        )

        # 4) Beads-subtracted normalized curves
        plot_sawtooth_grid_for_dye(
            tidy=tidy_bg,
            dye_name=dye_name,
            use_norm=True,
            out_dir=SAWTOOTH_GRID_DIR_BG_NORM,
            sample_order=sample_order,
        )

    print(f"\n   ✓ Raw grids written to:         {SAWTOOTH_GRID_DIR}")
    print(f"   ✓ Normalized grids written to:  {SAWTOOTH_GRID_DIR_NORM}")
    print(f"   ✓ BG-sub raw grids to:          {SAWTOOTH_GRID_DIR_BG}")
    print(f"   ✓ BG-sub norm grids to:         {SAWTOOTH_GRID_DIR_BG_NORM}")



if __name__ == "__main__":
    main()



# # ---------------------------------------------------------
# # CONFIGURATION
# # ---------------------------------------------------------

# DATA_DIR = "."  # folder containing your *ex.xlsx files
# WELL_COLUMN = "Well"
# WAVELENGTH_COLUMN = "Unnamed: 1"  # column that holds "Wavelength [nm]"
# RAW_FLAG = "Raw Data  (Em Spectrum)"  # value in WELL_COLUMN for real spectra

# # 24 samples: A01–F04
# SAMPLE_COLUMNS = [f"{c}{i:02d}" for c in "ABCDEF" for i in range(1, 5)]

# # # Metadata CSV mapping wells to amyloid type
# # METADATA_CSV = "sample_metadata.csv"
# DYES = ["NIAD4", "LDS698", "Bf188", "DQTCI", "MCAAD3", "ThT"]
# SAMPLES = ["Beads", "Monomeric Tau", "Tau fibrils", "a-synuclein fibrils"]

# RANDOM_STATE = 0

# # ---------------------------------------------------------
# # LOADING & BUILDING EMBER-LIKE PROFILES
# # ---------------------------------------------------------

# def load_excel_files(data_dir):
#     """Find all *.xlsx files and return as an ordered dict keyed by excitation name."""
#     paths = sorted(glob.glob(os.path.join(data_dir, "*.xlsx")))
#     ex_files = OrderedDict()
#     for path in paths:
#         # Use filename (without extension) as excitation label, e.g. "350ex"
#         label = os.path.splitext(os.path.basename(path))[0]
#         ex_files[label] = path
#     if not ex_files:
#         raise FileNotFoundError("No *ex.xlsx files found in DATA_DIR.")
#     return ex_files


# def extract_profiles_from_file(path):
#     """
#     From one excitation Excel file, return:
#       wavelengths: np.array of emission wavelengths
#       intensity_dict: {sample_id -> np.array(intensities)}
#     Assumes:
#       - rows with WELL_COLUMN == RAW_FLAG are actual spectra
#       - WAVELENGTH_COLUMN holds emission wavelength
#       - SAMPLE_COLUMNS hold intensity values
#     """
#     df = pd.read_excel(path)

#     # Keep only spectral rows
#     df_data = df[df[WELL_COLUMN] == RAW_FLAG].copy()

#     wavelengths = df_data[WAVELENGTH_COLUMN].astype(float).to_numpy()

#     intensity_dict = {}
#     for col in SAMPLE_COLUMNS:
#         intensity_dict[col] = df_data[col].astype(float).to_numpy()

#     return wavelengths, intensity_dict


# def build_ember_matrix(ex_files):
#     """
#     Combine all excitation files into:
#       X: (n_samples, n_features) EMBER-style profiles
#       sample_ids: list of sample IDs corresponding to rows of X
#       em_waves_per_ex: list of wavelength arrays (one per excitation)
#       ex_labels: list of excitation labels in the same order
#     """
#     sample_ids = SAMPLE_COLUMNS[:]  # 24 wells
#     per_sample_blocks = {sid: [] for sid in sample_ids}
#     em_waves_per_ex = []
#     ex_labels = []

#     for ex_label, path in ex_files.items():
#         wavelengths, intensities = extract_profiles_from_file(path)
#         em_waves_per_ex.append(wavelengths)
#         ex_labels.append(ex_label)

#         for sid in sample_ids:
#             per_sample_blocks[sid].append(intensities[sid])

#     # Concatenate excitation blocks for each sample
#     profiles = []
#     for sid in sample_ids:
#         profiles.append(np.concatenate(per_sample_blocks[sid]))
#     X = np.vstack(profiles)

#     # Min–max normalize each profile to [0, 1] (matching Yang et al.)
#     X_min = X.min(axis=1, keepdims=True)
#     X_max = X.max(axis=1, keepdims=True)
#     X_ember = (X - X_min) / (X_max - X_min + 1e-9)

#     return X_ember, sample_ids, em_waves_per_ex, ex_labels


# # ---------------------------------------------------------
# # PLOTTING: SAWTOOTH EMBER PROFILE
# # ---------------------------------------------------------

# def plot_sawtooth_profile(sample_index, X, em_waves_per_ex, ex_labels, title=None):
#     """
#     Plot the concatenated "sawtooth" profile for one sample.
#     """
#     profile = X[sample_index]

#     # Split the concatenated vector back into blocks
#     block_lengths = [len(w) for w in em_waves_per_ex]
#     blocks = []
#     start = 0
#     for L in block_lengths:
#         blocks.append(profile[start:start+L])
#         start += L

#     plt.figure(figsize=(10, 4))
#     offset = 0.0

#     for block_y, waves, ex_label in zip(blocks, em_waves_per_ex, ex_labels):
#         # Shift each block along x so they don't overlap – makes a sawtooth
#         x = waves + offset
#         plt.plot(x, block_y, marker=".", linewidth=1)
#         # Add a little gap between excitation blocks
#         offset = x[-1] + (waves[1] - waves[0]) * 2

#     plt.xlabel("Concatenated emission windows (nm, offset by excitation)")
#     plt.ylabel("Normalized intensity")
#     if title:
#         plt.title(title)
#     plt.tight_layout()
#     plt.show()


# # ---------------------------------------------------------
# # PCA, UMAP, & DISCRIMINATION SCORE
# # ---------------------------------------------------------

# def run_pca_umap(X):
#     """
#     Standardize features, run PCA & UMAP.
#     Returns:
#       pca_scores: (n_samples, 2)
#       umap_embedding: (n_samples, 2)
#     """
#     scaler = StandardScaler()
#     X_std = scaler.fit_transform(X)

#     pca = PCA(n_components=2, random_state=RANDOM_STATE)
#     pca_scores = pca.fit_transform(X_std)

#     reducer = umap.UMAP(random_state=RANDOM_STATE)
#     umap_embedding = reducer.fit_transform(X_std)

#     return pca_scores, umap_embedding


# def plot_embedding(embedding, labels, title, ax=None):
#     """
#     Scatter plot of PCA or UMAP embedding colored by amyloid type.
#     """
#     if ax is None:
#         fig, ax = plt.subplots(figsize=(6, 5))

#     labels = np.array(labels)
#     unique = sorted(np.unique(labels))

#     for lab in unique:
#         mask = labels == lab
#         ax.scatter(
#             embedding[mask, 0],
#             embedding[mask, 1],
#             label=lab,
#             s=60,
#             alpha=0.8
#         )

#     ax.set_xlabel(f"{title} 1")
#     ax.set_ylabel(f"{title} 2")
#     ax.set_title(title)
#     ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")


# def discrimination_score(features, labels, n_splits=5, n_repeats=10):
#     """
#     Approximate Yang et al.'s quadratic discrimination score:
#     train a Quadratic Discriminant Analysis classifier on the
#     2D embedding (PCA or UMAP) with repeated stratified CV
#     and report mean accuracy (%).:contentReference[oaicite:1]{index=1}
#     """
#     labels = np.array(labels)
#     qda = QuadraticDiscriminantAnalysis(store_covariance=True)

#     rng = np.random.default_rng(RANDOM_STATE)
#     all_scores = []

#     for _ in range(n_repeats):
#         cv = StratifiedKFold(
#             n_splits=n_splits,
#             shuffle=True,
#             random_state=int(rng.integers(0, 1_000_000))
#         )
#         scores = cross_val_score(qda, features, labels, cv=cv)
#         all_scores.extend(scores)

#     return 100 * np.mean(all_scores)


# # ---------------------------------------------------------
# # MAIN
# # ---------------------------------------------------------

# def main():
#     # 1) Load excitation files and build EMBER matrix
#     ex_files = load_excel_files(DATA_DIR)
#     print("Found excitation files:", list(ex_files.keys()))

#     X_ember, sample_ids, em_waves_per_ex, ex_labels = build_ember_matrix(ex_files)
#     print("EMBER matrix shape:", X_ember.shape)

#     # 2) Load amyloid labels
#     meta = pd.read_csv(METADATA_CSV, comment="#")
#     meta = meta.set_index("sample_id")["amyloid"].to_dict()

#     labels = [meta[sid] for sid in sample_ids]

#     # 3) Plot one example sawtooth profile (change index if you like)
#     plot_sawtooth_profile(
#         sample_index=0,
#         X=X_ember,
#         em_waves_per_ex=em_waves_per_ex,
#         ex_labels=ex_labels,
#         title=f"Sawtooth profile: {sample_ids[0]} ({labels[0]})",
#     )

#     # 4) PCA & UMAP
#     pca_scores, umap_embedding = run_pca_umap(X_ember)

#     fig, axes = plt.subplots(1, 2, figsize=(12, 5))
#     plot_embedding(pca_scores, labels, "PCA", ax=axes[0])
#     plot_embedding(umap_embedding, labels, "UMAP", ax=axes[1])
#     plt.tight_layout()
#     plt.show()

#     # 5) Discrimination scores
#     pca_disc = discrimination_score(pca_scores, labels)
#     umap_disc = discrimination_score(umap_embedding, labels)

#     print(f"PCA discrimination score:  {pca_disc:5.1f}%")
#     print(f"UMAP discrimination score: {umap_disc:5.1f}%")


# if __name__ == "__main__":
#     main()
