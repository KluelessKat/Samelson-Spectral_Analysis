#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
spectral_combine_csv.py

Combine sample rows (rows whose FIRST column looks like [amyloid-dye],
with no underscores) from multiple Spillover CSV files into one CSV.

Usage:
    # Option 1: use paths hard-coded in main()
    python spectral_combine_csv.py

    # Option 2: pass folders and/or CSV files from terminal
    python spectral_combine_csv.py /path/to/run1 /path/to/run2 /path/to/file3.csv
"""

import os
import sys
import re
import pandas as pd
import fnmatch

# def find_csv_files(folder_dirs=[], filename_pattern="*_A.csv"):
#     """
#     Iteratively walk root_dir and return full paths of all CSV files
#     whose names match filename_pattern (e.g. 'Spillover_LSM_A.csv').
#     """
#     csv_paths = []
#     if folder_dirs != None:
#         for folder in folder_dirs:
#             for dirpath, dirnames, filenames in os.walk(folder):
#                 for fname in filenames:
#                     if fnmatch.fnmatch(fname, filename_pattern):
#                         full_path = os.path.join(dirpath, fname)
#                         csv_paths.append(full_path)
#     else:
#         print("no folder directories provided")
#     return csv_paths


# -------------------------------------------------------------------
# Sample-label detection: [amyloid-dye] with NO underscores
# -------------------------------------------------------------------

SAMPLE_PATTERN = re.compile(r"^\[[^_\]]+-[^_\]]+\]$")

def is_sample_label(value: str) -> bool:
    """
    Return True if value looks like [amyloid-dye] with no underscores.
    Examples that match: [TauF-bf188], [asyn-ThT]
    Examples that do NOT match: [blank], [TauF_bf188], TauF-bf188
    """
    if value is None:
        return False
    s = str(value).strip()
    return bool(SAMPLE_PATTERN.match(s))

# -------------------------------------------------------------------
# Path handling
# -------------------------------------------------------------------

def expand_paths_to_csvs(paths):
    """
    For each input path:
      - If it's a directory → walk it recursively and collect all *.csv
      - If it's a file ending in .csv → include it
      - Otherwise → print a warning and skip

    Returns: sorted list of unique CSV file paths.
    """
    csv_files = []

    for p in paths:
        p = os.path.abspath(p)
        if os.path.isdir(p):
            # walk recursively through the folder
            for dirpath, dirnames, filenames in os.walk(p):
                for fname in filenames:
                    if fname.lower().endswith(".csv"):
                        csv_files.append(os.path.join(dirpath, fname))
        elif os.path.isfile(p) and p.lower().endswith(".csv"):
            csv_files.append(p)
        else:
            print(f"⚠️ Skipping non-CSV path: {p}")

    # deduplicate and sort
    csv_files = sorted(set(csv_files))
    return csv_files

# -------------------------------------------------------------------
# Core combiner
# -------------------------------------------------------------------

def combine_sample_rows(csv_paths, output_csv="combined_samples.csv"):
    """
    For each CSV in csv_paths:
      - read it
      - detect sample rows by inspecting the FIRST column
      - keep only rows whose first-column value looks like [amyloid-dye]
      - add SourceFile and SourceFolder columns
    Concatenate all sample rows and write to output_csv.
    """
    if not csv_paths:
        raise ValueError("No CSV files provided to combine.")

    dfs = []

    for path in csv_paths:
        print(f"📄 Reading: {path}")
        df = pd.read_csv(path, header=0)

        # The first column holds the row labels (Channels, Mask, [TauF-bf188], etc.)
        first_col = df.columns[0]
        labels = df[first_col].astype(str)

        # mask out only real sample rows: [amyloid-dye] with no underscores
        mask = labels.apply(is_sample_label)
        df_samples = df[mask].copy()

        if df_samples.empty:
            print(f"  ⚠️ No matching [amyloid-dye] rows found in this file.")
            continue

        # (optional) add a normalized Sample column for clarity
        df_samples.insert(1, "Sample", labels[mask].values)

        # provenance
        df_samples["SourceFile"] = os.path.basename(path)
        df_samples["SourceFolder"] = os.path.dirname(os.path.abspath(path))

        dfs.append(df_samples)

    if not dfs:
        raise RuntimeError("No valid sample rows found in any CSV.")

    combined = pd.concat(dfs, ignore_index=True)

    combined.to_csv(output_csv, index=False)
    print(f"\n✅ Combined samples written to: {output_csv}")
    print(f"   Rows: {len(combined)} from {len(csv_paths)} CSV files")

    return combined

# -------------------------------------------------------------------
# main()
# -------------------------------------------------------------------

def main():
    # If paths are given on the command line, use them
    if len(sys.argv) > 1:
        input_paths = sys.argv[1:]
    else:
        # Otherwise, manually specify folders and/or CSV files here:
        input_paths = [
            # Example: folder containing Spillover_LSM_A.csv etc.
            # "/Users/you/Downloads/Spectral_Data_Processing/a-syn,bf188,tht_20251115180949",
            # Example: explicit CSV files
            "/Users/katherinezhang/Downloads/Spectral_Data_Processing/a-syn,bf188,tht_20251115180949/Spillover_LSM_A.csv",
            "/Users/katherinezhang/Downloads/Spectral_Data_Processing/mTau,Tau-bf188,ThT_20251119170806/Spillover_LSM_A.csv",
        ]

    csv_paths = expand_paths_to_csvs(input_paths)

    if not csv_paths:
        print("❌ No CSV files found. Edit input_paths in main() or pass paths via terminal.")
        return

    print("Found CSV files:")
    for p in csv_paths:
        print("  ", p)

    combine_sample_rows(csv_paths, output_csv="combined_samples.csv")


if __name__ == "__main__":
    main()
