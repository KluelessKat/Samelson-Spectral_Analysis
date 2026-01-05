# -*- coding: utf-8 -*-
"""
spectral_preprocessing.py

Input:  combined_samples.csv (from spectral_combine_csv.py)
Output:
  - tidy Excel file with columns:
      Sample, Amyloid, Dye, Excitation_nm,
      Detector_CH, Emission_nm_start, Intensity, OrderInBlock
  - per-dye sawtooth PNGs:
      one figure per dye, rows = amyloids, cols = excitations

        Sawtooth plot directory format example:
        combined_samples_sawtooth_plots_20251222_143501/
        ├── uncollapsed_raw/
        ├── uncollapsed_norm/
        ├── collapsed_raw/
        ├── collapsed_globalnorm_stacked/
        └── collapsed_staggered_norm/
"""


import os
from io import BytesIO
import xlsxwriter
from glob import glob
import fnmatch
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors


"""Instrument-specific configuration"""

# SONY Channel start wavelengths (nm) for CH1–CH35
CH_START_NM = {
    1: 360.5,
    2: 418.5,
    3: 446,
    4: 493.9,
    5: 506.1,
    6: 518.3,
    7: 530.3,
    8: 542.3,
    9: 554.2,
    10: 566.1,
    11: 577.8,
    12: 589.5,
    13: 601.1,
    14: 612.6,
    15: 624,
    16: 635.3,
    17: 646.6,
    18: 657.7,
    19: 668.8,
    20: 679.8,
    21: 690.8,
    22: 701.6,
    23: 712.4,
    24: 723,
    25: 733.6,
    26: 744.2,
    27: 754.6,
    28: 765,
    29: 775.2,
    30: 785.4,
    31: 795.5,
    32: 805.5,
    33: 815.5,
    34: 825.3,
    35: 835.1,
}

# Block layout (channel index ranges) – 1-based indices as in your CSV header row
BLOCKS = [
    (1, 35),
    (38, 72),
    (75, 109),
    (115, 146),
    (158, 183),
    (202, 220),
]

# Excitation wavelengths, in the same order as BLOCKS
EXCITATIONS = [320, 355, 405, 488, 561, 637]

# Metadata labels to ignore as “samples”
BAD_LABELS = {
    "Channels",
    "Mask",
    "Offset",
    "PDF_Std",
    "Noise_p",
    "AutoFluorescence_m",
    "AutoFluorescence",
    "AutoFluorescence_p",
}

# -------------------------------------------------------------------
# HELPER FUNCTIONS
# -------------------------------------------------------------------

def drop_leading_zeros(seg: pd.DataFrame, eps: float = 1e-6) -> pd.DataFrame:
    """Return seg with leading ~zero intensity points removed (for plotting only)."""
    vals = seg["Intensity"].to_numpy(dtype=float)
    mask = np.abs(vals) > eps
    if not mask.any():
        return seg.iloc[0:0].copy()  # all zeros → empty
    first = int(np.argmax(mask))
    return seg.iloc[first:].copy()

def safe_name(text: str) -> str:
    """Make a string safe for filenames."""
    bad = '[]:*?/\\'
    s = "".join(ch for ch in str(text) if ch not in bad)
    s = s.replace(" ", "_")
    return s

def make_plot_dirs(base: str):
    """
    Create a timestamped plot root folder plus subfolders for each plot type.
    Returns a dict of {key: folder_path}.
    """
    #ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    ts = ""
    root = f"{base}_sawtooth_plots_{ts}"

    subdirs = {
        "uncollapsed_raw": os.path.join(root, "uncollapsed_raw"),
        "uncollapsed_norm": os.path.join(root, "uncollapsed_norm"),
        "collapsed_raw": os.path.join(root, "collapsed_raw"),
        "collapsed_globalnorm": os.path.join(root, "collapsed_globalnorm_stacked"),
        "collapsed_staggered": os.path.join(root, "collapsed_staggered_norm"),
    }

    for p in subdirs.values():
        os.makedirs(p, exist_ok=True)

    return root, subdirs

def build_tidy_from_combined_csv(csv_path: str) -> pd.DataFrame:
    """
    Read combined_samples.csv and build a tidy spectral dataframe.

    We:
      - use the 'Sample' column as label
      - treat only columns whose names are integers as detector channels
      - map channels into Detector_CH and Emission_nm_start using BLOCKS/CH_START_NM
    """
    df = pd.read_csv(csv_path, header=0)

    if "Sample" not in df.columns:
        raise ValueError("combined CSV must have a 'Sample' column.")

    label_col = "Sample"

    # numeric-channel columns only (e.g. '1', '2', ..., '259')
    channel_cols = []
    header_map = {}
    for c in df.columns:
        if c == label_col:
            continue
        if c in ["SourceFile", "SourceFolder", "0"]:
            continue
        try:
            ci = int(str(c).strip())
            header_map[ci] = c
            channel_cols.append(c)
        except Exception:
            # metadata columns like Unnamed:260 etc.
            continue

    if not header_map:
        raise RuntimeError("No numeric channel columns found in combined CSV.")

    tidy_records = []

    for _, row in df.iterrows():
        sample_name = str(row[label_col]).strip()

        for block_idx, (start, end) in enumerate(BLOCKS):
            ex = EXCITATIONS[block_idx]
            block_channels = list(range(start, end + 1))
            width = len(block_channels)
            start_ch = 35 - width + 1
            ch_indices = list(range(start_ch, 35 + 1))

            for ch_idx_in_block, chan_idx in zip(ch_indices, block_channels):
                colname = header_map.get(chan_idx)
                if colname is None:
                    val = np.nan
                else:
                    v = row.get(colname, np.nan)
                    try:
                        val = float(str(v).replace(",", ""))
                    except Exception:
                        val = np.nan

                wl = CH_START_NM.get(ch_idx_in_block, np.nan)

                tidy_records.append(
                    {
                        "Sample": sample_name,
                        "Excitation_nm": ex,
                        "ChannelIndex": chan_idx,
                        "Detector_CH": ch_idx_in_block,
                        "Emission_nm_start": wl,
                        "Intensity": val,
                    }
                )

    tidy = pd.DataFrame(tidy_records)
    if tidy.empty:
        raise RuntimeError("No records generated from combined CSV.")

    tidy["OrderInBlock"] = tidy.groupby(["Sample", "Excitation_nm"]).cumcount() + 1

    # parse Amyloid / Dye from Sample: assume "[Amyloid-Dye]"
    tidy["SampleClean"] = tidy["Sample"].str.strip("[]")
    tidy[["Amyloid", "Dye"]] = (tidy["SampleClean"].str.split("-", n=1, expand=True)
    .apply(lambda col: col.str.strip().str.title())
    )

    return tidy

def add_normalized_intensity(tidy: pd.DataFrame) -> pd.DataFrame:
    """
    Add Intensity_norm column:
      for each (Amyloid, Dye, Excitation_nm) group,
      divide intensities by the max in that group.
    """
    tidy = tidy.copy()
    def _norm(x):
        m = x.max()
        return x / m if m > 0 else 0.0

    tidy["Intensity_norm"] = (
        tidy.groupby(["Amyloid", "Dye", "Excitation_nm"])["Intensity"]
            .transform(_norm)
    )
    return tidy

# -------------------------------------------------------------------
# PLOT
# -------------------------------------------------------------------
def plot_sawtooths_per_dye(tidy: pd.DataFrame, dye_name: str, outdir: str,
                           normalize: bool = False):
    """
    Plot per-dye sawtooths:
      rows = amyloids, cols = excitations
      option to normalize each excitation's curve to max = 1
    """
    sub = tidy[tidy["Dye"] == dye_name].copy()
    if sub.empty:
        print(f"No rows for dye {dye_name}")
        return

    amyloids = sorted(sub["Amyloid"].unique())

    # Use only excitations that actually exist for this dye (safer than fixed EXCITATIONS)
    excitations_here = [ex for ex in EXCITATIONS if ex in set(sub["Excitation_nm"].dropna().unique())]
    if not excitations_here:
        print(f"No excitations found for dye {dye_name}")
        return
    
    n_rows = len(amyloids)
    n_cols = len(excitations_here)

    cmap = cm.get_cmap("viridis", len(excitations_here))
    ex_to_color = {ex: cmap(i) for i, ex in enumerate(excitations_here)}

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.2 * n_cols, 3.0 * n_rows),
        sharex=True, sharey=True
    )

    if n_rows == 1:
        axes = np.expand_dims(axes, axis=0)  # shape (1, n_cols)

    for i, amy in enumerate(amyloids):
        sub_amy = sub[sub["Amyloid"] == amy]
        for j, ex in enumerate(excitations_here):
            ax = axes[i, j]
            seg = sub_amy[sub_amy["Excitation_nm"] == ex].copy()
            seg = seg.sort_values("Emission_nm_start")
            seg = drop_leading_zeros(seg)
            if seg.empty:
                ax.set_visible(False)
                continue

            x = seg["Emission_nm_start"].to_numpy(dtype=float)
            y = seg["Intensity"].to_numpy(dtype=float)
            if normalize and y.max() > 0:
                y = y / y.max()

            ax.plot(x, y, marker="o", linewidth=2, color=ex_to_color[ex],)

            if i == 0:
                ax.set_title(f"{ex} nm")
            if j == 0:
                ylabel = amy
                if normalize:
                    ylabel += "\n(norm. max = 1)"
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(axis="y", left=False, labelleft=False)

    title = f"Uncollapsed Sawtooths – {dye_name}"
    if normalize:
        title += " (Normalized)"
    fig.suptitle(title, fontsize=16, y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    #Save Plots
    os.makedirs(outdir, exist_ok=True)
    fname = f"uncollapsed_sawtooths_{dye_name}{'_norm' if normalize else ''}.png"
    save_path = os.path.join(outdir, fname)
    plt.savefig(save_path, dpi=300)
    plt.close(fig)
    print("Saved plot:", save_path)

def plot_collapsed_spectra_raw(tidy: pd.DataFrame, dye_name: str, outdir: str):
    """
    For a given dye, plot a collapsed spectrum for each amyloid:
      - x-axis: Emission_nm_start
      - y-axis: raw Intensity
      - one panel = one (Amyloid, Dye) combination
      - multiple excitation curves overlaid, like the EMBER example.

    We DO NOT pad missing emission wavelengths; we just plot the
    wavelengths we actually have (matplotlib ignores NaNs).
    """
    sub = tidy[tidy["Dye"] == dye_name].copy()
    if sub.empty:
        print(f"No rows for dye {dye_name} (raw collapsed).")
        return

    # use only excitations that actually appear for this dye
    excitations_here = sorted(sub["Excitation_nm"].dropna().unique())
    cmap = cm.get_cmap("viridis", len(excitations_here))

    amyloids = sorted(sub["Amyloid"].dropna().unique())
    if not amyloids:
        print(f"No amyloids found for dye {dye_name} (raw collapsed).")
        return

    os.makedirs(outdir, exist_ok=True)

    for amy in amyloids:
        sub_amy = sub[sub["Amyloid"] == amy].copy()
        if sub_amy.empty:
            continue

        fig, ax = plt.subplots(figsize=(5, 4))

        for i, ex in enumerate(excitations_here):
            seg = sub_amy[sub_amy["Excitation_nm"] == ex].copy()
            if seg.empty:
                continue

            seg = seg.sort_values("Emission_nm_start")
            seg = drop_leading_zeros(seg)  # remove leading ~0’s only
            if seg.empty:
                continue

            x = seg["Emission_nm_start"].to_numpy(dtype=float)
            y = seg["Intensity"].to_numpy(dtype=float)

            # No padding; just plot existing x,y
            ax.plot(x, y, marker="o", linewidth=2, color=cmap(i), label=f"{ex} nm")

        ax.set_xlabel("Emission (nm)")
        ax.set_ylabel("Intensity (A.U.)")
        ax.set_title(f"{amy} – {dye_name} (collapsed, raw)")
        ax.legend(title="Excitation nm", fontsize=8)
        ax.grid(False)

        plt.tight_layout()

        fname = f"collapsed_raw_{safe_name(dye_name)}_{safe_name(amy)}.png"
        save_path = os.path.join(outdir, fname)
        plt.savefig(save_path, dpi=300)
        plt.close(fig)
        print("Saved plot:", save_path)

def plot_collapsed_spectra_global_norm(
    tidy: pd.DataFrame,
    dye_name: str,
    outdir: str,
    offset_step: float = 0.5,
):
    """
    For a given dye and each amyloid, plot collapsed spectra where:

      - each excitation curve is normalized to ITS OWN max
        (so each has a peak of 1.0 before shifting),
      - curves are then vertically offset by `offset_step` increments:
          0*offset_step for the first excitation,
          1*offset_step for the second,
          2*offset_step for the third, ...

    This mimics the stacked presentation in the example figure:
      x-axis: Emission_nm_start
      y-axis: normalized intensity + offset
    """
    sub = tidy[tidy["Dye"] == dye_name].copy()
    if sub.empty:
        print(f"No rows for dye {dye_name} (per-curve norm + offset).")
        return

    excitations_here = sorted(sub["Excitation_nm"].dropna().unique())
    if not excitations_here:
        print(f"No excitations found for dye {dye_name} (per-curve norm + offset).")
        return

    k = len(excitations_here)
    cmap = cm.get_cmap("viridis", k)

    amyloids = sorted(sub["Amyloid"].dropna().unique())
    if not amyloids:
        print(f"No amyloids found for dye {dye_name} (per-curve norm + offset).")
        return

    os.makedirs(outdir, exist_ok=True)

    for amy in amyloids:
        sub_amy = sub[sub["Amyloid"] == amy].copy()
        if sub_amy.empty:
            continue

        fig, ax = plt.subplots(figsize=(5, 4))

        max_y_seen = 0.0

        for i, ex in enumerate(excitations_here):
            seg = sub_amy[sub_amy["Excitation_nm"] == ex].copy()
            if seg.empty:
                continue

            seg = seg.sort_values("Emission_nm_start")
            seg = drop_leading_zeros(seg)
            if seg.empty:
                continue

            x = seg["Emission_nm_start"].to_numpy(dtype=float)
            y_raw = seg["Intensity"].to_numpy(dtype=float)

            curve_max = np.max(y_raw)
            if not np.isfinite(curve_max) or curve_max <= 0:
                # skip excitation with no positive signal
                continue

            # normalize this excitation to its own max → peak at 1.0
            y_norm = y_raw / curve_max

            # vertical offset to “stack” the curves
            offset = i * offset_step
            y = y_norm + offset

            max_y_seen = max(max_y_seen, y.max())

            ax.plot(
                x,
                y,
                marker="o",
                linewidth=2,
                color=cmap(i),
                label=f"{ex} nm",
            )

        if max_y_seen == 0:
            # nothing actually plotted
            plt.close(fig)
            continue

        ax.set_xlabel("Emission (nm)")
        ax.set_ylabel("Normalized intensity (per excitation) + offset")
        ax.set_ylim(0, max_y_seen + 0.1)

        ax.set_title(f"{amy} – {dye_name} (collapsed, each curve max=1, stacked)")
        ax.legend(title="Excitation nm", fontsize=8)
        ax.grid(False)

        plt.tight_layout()
        fname = f"collapsed_globalnorm_{safe_name(dye_name)}_{safe_name(amy)}.png"
        save_path = os.path.join(outdir, fname)
        plt.savefig(save_path, dpi=300)
        plt.close(fig)
        print("Saved plot:", save_path)


def plot_collapsed_spectra_staggered_norm(tidy: pd.DataFrame, dye_name: str, outdir: str):
    """
    For a given dye and amyloid, plot collapsed spectra where:

      - each excitation curve is normalized by ITS OWN max,
      - but scaled so that the curve maxima form a descending ladder
        starting near 1.0.

    Example:
      Ex list (sorted): [320, 355, 405, 488, 561, 637]
      k = 6
      step = 0.9 / (k - 1) -> 0.18
      target maxima: [1.0, 0.82, 0.64, 0.46, 0.28, 0.10]
    """
    sub = tidy[tidy["Dye"] == dye_name].copy()
    if sub.empty:
        print(f"No rows for dye {dye_name} (staggered-norm collapsed).")
        return

    excitations_here = sorted(sub["Excitation_nm"].dropna().unique())
    k = len(excitations_here)
    if k == 0:
        print(f"No excitations found for dye {dye_name} (staggered-norm).")
        return

    cmap = cm.get_cmap("viridis", k)

    # define target maxima descending from 1.0
    if k == 1:
        target_maxima = [1.0]
    else:
        step = 0.9 / (k - 1)   # from 1.0 down to ~0.1
        target_maxima = [1.0 - i * step for i in range(k)]

    amyloids = sorted(sub["Amyloid"].dropna().unique())
    if not amyloids:
        print(f"No amyloids found for dye {dye_name} (staggered-norm collapsed).")
        return

    os.makedirs(outdir, exist_ok=True)

    for amy in amyloids:
        sub_amy = sub[sub["Amyloid"] == amy].copy()
        if sub_amy.empty:
            continue

        fig, ax = plt.subplots(figsize=(5, 4))

        for i, ex in enumerate(excitations_here):
            seg = sub_amy[sub_amy["Excitation_nm"] == ex].copy()
            if seg.empty:
                continue

            seg = seg.sort_values("Emission_nm_start")
            seg = drop_leading_zeros(seg)
            if seg.empty:
                continue

            x = seg["Emission_nm_start"].to_numpy(dtype=float)
            y_raw = seg["Intensity"].to_numpy(dtype=float)
            curve_max = np.max(y_raw)

            if not np.isfinite(curve_max) or curve_max <= 0:
                # nothing to plot for this excitation
                continue

            desired_max = target_maxima[i]
            scale = desired_max / curve_max
            y = y_raw * scale

            ax.plot(x, y, marker="o", linewidth=2, color=cmap(i), label=f"{ex} nm (max={desired_max:.2f})")

        ax.set_xlabel("Emission (nm)")
        ax.set_ylabel("Staggered normalized intensity (A.U.)")
        ax.set_ylim(bottom=0)
        ax.set_title(f"{amy} – {dye_name} (collapsed, staggered norm)")
        ax.legend(title="Excitation nm", fontsize=7)
        ax.grid(False)

        plt.tight_layout()
        fname = f"collapsed_staggered_{safe_name(dye_name)}_{safe_name(amy)}.png"
        save_path = os.path.join(outdir, fname)
        plt.savefig(save_path, dpi=300)
        plt.close(fig)
        print("Saved plot:", save_path)

# -------------------------------------------------------------------
# HANDLE FILE EXPORTS
# -------------------------------------------------------------------

def export_tidy_excel(tidy: pd.DataFrame, excel_path: str):
    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        tidy.to_excel(writer, sheet_name="tidy", index=False)
    print("Wrote tidy Excel:", excel_path)

def export_sawtooths_by_dye_excel(tidy: pd.DataFrame, excel_path: str):
    """
    Create a single Excel file where:
      - each sheet corresponds to a dye
      - each sheet contains all rows for that dye, with:
          Sample, Amyloid, Dye, Excitation_nm,
          Emission_nm_start, Intensity, Intensity_norm, Detector_CH, OrderInBlock
      - rows are sorted by Amyloid, Excitation_nm, Emission_nm_start
    """
    dyes = sorted(tidy["Dye"].dropna().unique())
    if not dyes:
        print("No dyes found in tidy dataframe.")
        return

    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        for dye in dyes:
            sub = tidy[tidy["Dye"] == dye].copy()
            if sub.empty:
                continue

            sub = sub.sort_values(
                ["Amyloid", "Excitation_nm", "Emission_nm_start"]
            )

            # pick nice column order
            cols = [
                "Sample",
                "Amyloid",
                "Dye",
                "Excitation_nm",
                "Emission_nm_start",
                "Detector_CH",
                "OrderInBlock",
                "Intensity",
                "Intensity_norm",
            ]
            cols = [c for c in cols if c in sub.columns]
            sub = sub[cols]

            # Excel sheet names must be <=31 chars and no []:/\?*
            safe = "".join(ch for ch in str(dye) if ch not in '[]:*?/\\')
            sheet_name = safe[:31] or "Sheet"

            sub.to_excel(writer, sheet_name=sheet_name, index=False)

    print("Wrote sawtooth-by-dye Excel:", excel_path)


# -------- main --------

def main():
    # combined_samples.csv from script 1
    if len(sys.argv) > 1:
        csv_path = sys.argv[1]
    else:
        csv_path = "combined_samples.csv"   # edit if needed

    if not os.path.isfile(csv_path):
        print(f"❌ combined CSV not found: {csv_path}")
        return

    tidy = build_tidy_from_combined_csv(csv_path)
    #normalize intensities for each excitation-emission curve
    tidy_norm = add_normalized_intensity(tidy)

    print("Tidy shape:", tidy.shape)
    print("Dyes:", sorted(tidy["Dye"].dropna().unique()))
    print("Amyloids:", sorted(tidy["Amyloid"].dropna().unique()))

    base, _ = os.path.splitext(csv_path)
    
    # tidy_excel_path = base + "_tidy.xlsx"
    tidy_norm_excel_path = base + "_tidy_norm.xlsx"
    # export_tidy_excel(tidy, tidy_excel_path)
    export_tidy_excel(tidy_norm, tidy_norm_excel_path)
    
    # excel_saw_path = base + "_sawtooths_by_dye.xlsx"
    excel_saw_path_norm = base + "_sawtooths_by_dye_norm.xlsx"
    # export_sawtooths_by_dye_excel(tidy, excel_saw_path)
    export_sawtooths_by_dye_excel(tidy_norm, excel_saw_path_norm)

    dyes = sorted(tidy["Dye"].dropna().unique())

    plot_root, plot_dirs = make_plot_dirs(base)
    print(f"📁 Writing plots to: {plot_root}")

    for dye in dyes:
        # Uncollapsed sawtooths (rows=amyloids, cols=excitations)
        plot_sawtooths_per_dye(tidy, dye, plot_dirs["uncollapsed_raw"], normalize=False)
        plot_sawtooths_per_dye(tidy, dye, plot_dirs["uncollapsed_norm"], normalize=True)

        # Collapsed spectra per (amyloid, dye)
        plot_collapsed_spectra_raw(tidy, dye, plot_dirs["collapsed_raw"])
        plot_collapsed_spectra_global_norm(tidy, dye, plot_dirs["collapsed_globalnorm"])
        plot_collapsed_spectra_staggered_norm(tidy, dye, plot_dirs["collapsed_staggered"])


if __name__ == "__main__":
    main()


# def build_tidy_from_csv(csv_path: str) -> pd.DataFrame:
#     """Load an ID7000 Spillover_LSM_A-style CSV and return a tidy dataframe.

#     Columns of tidy:
#         Sample, Excitation_nm, ChannelIndex, Detector_CH,
#         Emission_nm_start, Intensity, OrderInBlock
#     """
#     df = pd.read_csv(csv_path, dtype=str, header=0)
#     first_col = df.columns[0]
#     df = df.rename(columns={first_col: "Label"})
#     num_cols = [c for c in df.columns if c != "Label"]

#     # Map 1-based channel index -> column name
#     header_map = {}
#     for i, c in enumerate(num_cols, start=1):
#         try:
#             ci = int(str(c).strip())
#             header_map[ci] = c
#         except Exception:
#             header_map[i] = c

#     tidy_records = []

#     for _, row in df.iterrows():
#         label = str(row["Label"]).strip()
#         if label in BAD_LABELS or label == "" or label.lower() == "nan":
#             continue

#         # drop _m / _p measurement rows
#         if label.endswith("_m") or label.endswith("_p") or label.endswith("_m]") or label.endswith("_p]"):
#             continue

#         sample_name = label

#         for block_idx, (start, end) in enumerate(BLOCKS):
#             ex = EXCITATIONS[block_idx]

#             # All channel indices (column positions) for this block in the CSV
#             block_channels = list(range(start, end + 1))
#             width = len(block_channels)

#             # Map these 'width' detectors onto the LAST 'width' spectral channels (CH1–35)
#             # e.g. width=32 → CH4–35, width=26 → CH10–35, width=19 → CH17–35
#             start_ch = 35 - width + 1
#             ch_indices = list(range(start_ch, 35 + 1))  # inclusive of 35

#             for ch_idx_in_block, chan_idx in zip(ch_indices, block_channels):
#                 colname = header_map.get(chan_idx)
#                 if colname is None:
#                     val = np.nan
#                 else:
#                     v = row.get(colname, np.nan)
#                     try:
#                         val = float(str(v).replace(",", ""))
#                     except Exception:
#                         val = np.nan

#                 wl = CH_START_NM.get(ch_idx_in_block, np.nan)

#                 tidy_records.append(
#                     {
#                         "Sample": sample_name,
#                         "Excitation_nm": ex,
#                         "ChannelIndex": chan_idx,
#                         "Detector_CH": ch_idx_in_block,
#                         "Emission_nm_start": wl,
#                         "Intensity": val,
#                     }
#                 )

#     tidy = pd.DataFrame(tidy_records)
#     if tidy.empty:
#         raise RuntimeError("No sample rows were parsed; check BAD_LABELS and CSV format.")

#     tidy["OrderInBlock"] = tidy.groupby(["Sample", "Excitation_nm"]).cumcount() + 1
#     return tidy


# def build_sawtooth_uncollapsed_df(tidy: pd.DataFrame, sample_name: str) -> pd.DataFrame:
#     """Build uncollapsed saw-tooth style data for one sample."""
#     sub = tidy[tidy["Sample"] == sample_name].copy()
#     dfs = []
#     offset = 0
#     for ex in EXCITATIONS:
#         seg = sub[sub["Excitation_nm"] == ex].copy()
#         seg = seg.sort_values("Emission_nm_start")
#         seg = drop_leading_zeros(seg)
#         if seg.empty:
#             continue
#         seg["SawX"] = np.arange(1, len(seg) + 1) + offset
#         seg["SegmentLabel"] = f"{ex}nm"
#         offset += len(seg) + 3
#         dfs.append(
#             seg[
#                 [
#                     "SawX",
#                     "Emission_nm_start",
#                     "Intensity",
#                     "Excitation_nm",
#                     "SegmentLabel",
#                     "OrderInBlock",
#                 ]
#             ]
#         )
#     return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

# def build_per_excitation_dfs(tidy: pd.DataFrame, sample_name: str, drop_zeros: bool = True):
#     """
#     For a given sample, return a dict: {excitation_nm: DataFrame}
#     Each DataFrame has columns: Emission_nm_start, Intensity, Detector_CH, Excitation_nm.
#     Optionally drops leading ~zero-intensity points.
#     """
#     result = {}
#     sub = tidy[tidy["Sample"] == sample_name].copy()

#     for ex in EXCITATIONS:
#         seg = sub[sub["Excitation_nm"] == ex].copy()
#         seg = seg.sort_values("Emission_nm_start")
#         if drop_zeros:
#             seg = drop_leading_zeros(seg)
#         if seg.empty:
#             continue
#         # Keep only the columns we usually care about
#         result[ex] = seg[["Emission_nm_start", "Intensity", "Detector_CH", "Excitation_nm"]]

#     return result

# """## Load CSV and build tidy spectral dataframe
# Edit `csv_path` below to point to your `Spillover_LSM_A.csv` file, then run the cell.
# """

# # ---------- Plotting ----------

# def plot_side_by_side_per_excitation(tidy: pd.DataFrame, sample_name: str, outdir: str):
#     sub = tidy[tidy["Sample"] == sample_name].copy()
#     if sub.empty:
#         return

#     # one panel per excitation
#     n_exc = len(EXCITATIONS)
#     fig, axes = plt.subplots(1, n_exc, figsize=(3 * n_exc, 4), sharey=True)

#     if n_exc == 1:
#         axes = [axes]

#     # global y-limits
#     ymin, ymax = np.inf, -np.inf
#     per_exc = {}
#     for ex in EXCITATIONS:
#         seg = sub[sub["Excitation_nm"] == ex].copy()
#         seg = seg.sort_values("Emission_nm_start")
#         seg = drop_leading_zeros(seg)
#         if seg.empty:
#             continue
#         per_exc[ex] = seg
#         ymin = min(ymin, seg["Intensity"].min())
#         ymax = max(ymax, seg["Intensity"].max())

#     if not per_exc:
#         return

#     pad = 0.05 * (ymax - ymin) if ymax > ymin else 1.0
#     ymin -= pad
#     ymax += pad

#     for i, ex in enumerate(EXCITATIONS):
#         ax = axes[i]
#         seg = per_exc.get(ex)
#         if seg is not None:
#             ax.plot(seg["Emission_nm_start"], seg["Intensity"], linewidth=2)
#         ax.set_title(f"{ex} nm")
#         ax.set_xlabel("λ (nm)")
#         if i == 0:
#             ax.set_ylabel("Intensity (A.U.)")
#         else:
#             ax.tick_params(axis="y", left=False, labelleft=False)
#         ax.set_ylim(ymin, ymax)

#     plt.suptitle(f"{sample_name} – Per-Excitation Spectra", y=1.02)
#     plt.tight_layout()

#     os.makedirs(outdir, exist_ok=True)
#     safe = "".join(ch for ch in sample_name if ch not in '[]:*?/\\')
#     png_path = os.path.join(outdir, f"{safe}_per_excitation.png")
#     plt.savefig(png_path, dpi=300)
#     plt.close(fig)
#     print("Saved plot:", png_path)

# # ---------- Export ----------

# def export_tidy_to_excel(tidy: pd.DataFrame, excel_path: str):
#     with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
#         tidy.to_excel(writer, sheet_name="tidy", index=False)
#     print("Wrote tidy Excel:", excel_path)

# # ---------- Main ----------

# def main():
#     # 1) Get CSV path from terminal or use default
#     if len(sys.argv) > 1:
#         csv_path = sys.argv[1]
#     else:
#         csv_path = "Spillover_LSM_A.csv"  # <-- edit this

#     if not os.path.isfile(csv_path):
#         print(f"❌ CSV not found: {csv_path}")
#         return

#     # 2) Build tidy dataframe
#     tidy = build_tidy_from_csv(csv_path)
#     samples = sorted(tidy["Sample"].unique())
#     print(f"Loaded tidy dataframe with shape: {tidy.shape}")
#     print("Samples:", samples)

#     # 3) Export tidy Excel
#     base, _ = os.path.splitext(csv_path)
#     excel_path = base + "_tidy.xlsx"
#     export_tidy_to_excel(tidy, excel_path)

#     # 4) Generate and save side-by-side spectra plots
#     outdir = base + "_plots"
#     for sample in samples:
#         plot_side_by_side_per_excitation(tidy, sample, outdir)


# if __name__ == "__main__":
#     main()