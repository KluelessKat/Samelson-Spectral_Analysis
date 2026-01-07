#!/usr/bin/env python3
"""
aurora_rename_channels.py  (enhanced)

Outputs in --outdir:
1) Cleaned wide CSVs (per input):
   <stem>_cleaned_wide.csv
   - headers: channel token -> emission start (nm)
   - first data row inserted: excitation wavelength per channel-family

2) Combined tidy Excel:
   <name>_tidy_norm.xlsx   (sheet: tidy)

3) By-dye Excel:
   <name>_sawtooths_by_dye_norm.xlsx (one sheet per dye)

4) Sawtooth plots directory (if not --no-plots):
   <name>_sawtooth_plots_<timestamp>/
     ├── uncollapsed_raw/
     ├── uncollapsed_norm/
     ├── collapsed_raw/
     ├── collapsed_globalnorm_stacked/
     └── collapsed_staggered_norm/

Usage:
  python aurora_rename_channels.py --input-folder ./csvs --outdir ./out
  python aurora_rename_channels.py --inputs a.csv b.csv --outdir ./out --name run1
  python aurora_rename_channels.py --inputs a.csv --outdir ./out --unique-names
"""

from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Channel -> emission
# ---------------------------------------------------------------------

#"Wavelength Start (nm)"

CHANNEL_TO_START_NM: Dict[str, float] = {
    # Ultraviolet (UV)
    "UV1": 365, "UV2": 380, "UV3": 420, "UV4": 436, "UV5": 451, "UV6": 466,
    "UV7": 500, "UV8": 528, "UV9": 566, "UV10": 597, "UV11": 651, "UV12": 678,
    "UV13": 706, "UV14": 735, "UV15": 765, "UV16": 795,

    # Violet (V)
    "V1": 420, "V2": 436, "V3": 451, "V4": 466, "V5": 498, "V6": 516, "V7": 533,
    "V8": 571, "V9": 588, "V10": 605, "V11": 651, "V12": 678, "V13": 706,
    "V14": 735, "V15": 765, "V16": 795,

    # Blue (B)
    "B1": 498, "B2": 516, "B3": 533, "B4": 571, "B5": 588, "B6": 605,
    "B7": 653, "B8": 670, "B9": 688, "B10": 707, "B11": 728, "B12": 749,
    "B13": 772, "B14": 795,

    # Yellow-Green (YG)
    "YG1": 567, "YG2": 588, "YG3": 605, "YG4": 653, "YG5": 670, "YG6": 688,
    "YG7": 706, "YG8": 735, "YG9": 765, "YG10": 795,

    # Red (R)
    "R1": 653, "R2": 670, "R3": 688, "R4": 707, "R5": 728, "R6": 749,
    "R7": 772, "R8": 795,
}

#Center
CHANNEL_NM: Dict[str, float] = {
    # Ultraviolet (UV)
    "UV1": 373,  "UV2": 390,  "UV3": 429,  "UV4": 445,
    "UV5": 460,  "UV6": 475,  "UV7": 514,  "UV8": 542,
    "UV9": 581,  "UV10": 615, "UV11": 660, "UV12": 689,
    "UV13": 718, "UV14": 747, "UV15": 777, "UV16": 807,

    # Violet (V)
    "V1": 429,   "V2": 445,   "V3": 460,   "V4": 475,
    "V5": 508,   "V6": 525,   "V7": 542,   "V8": 581,
    "V9": 598,   "V10": 615,  "V11": 660,  "V12": 689,
    "V13": 718,  "V14": 747,  "V15": 777,  "V16": 807,

    # Blue (B)
    "B1": 508,   "B2": 525,   "B3": 542,   "B4": 581,
    "B5": 598,   "B6": 615,   "B7": 661,   "B8": 679,
    "B9": 697,   "B10": 717,  "B11": 738,  "B12": 760,
    "B13": 783,  "B14": 812,

    # Yellow-Green (YG)
    "YG1": 576,  "YG2": 598,  "YG3": 615,  "YG4": 661,
    "YG5": 679,  "YG6": 697,  "YG7": 718,  "YG8": 747,
    "YG9": 777,  "YG10": 807,

    # Red (R)
    "R1": 661,   "R2": 679,   "R3": 697,   "R4": 717,
    "R5": 738,   "R6": 760,   "R7": 783,   "R8": 812,
}


# Laser family -> excitation wavelength (nm)
LASER_EXCITATION_NM: Dict[str, int] = {
    "UV": 355,
    "V": 405,
    "B": 488,
    "YG": 561,
    "R": 640,
}

# Preferred excitation order for plotting
PREFERRED_EXCITATION_ORDER = [355, 405, 488, 561, 640]

# Regex to find channel tokens in column names
CHANNEL_TOKEN_RE = re.compile(
    r"\b(?P<prefix>UV|V|B|YG|R)(?P<num>\d{1,2})(?P<suffix>-[A-Za-z0-9]+)?\b"
)

BAD_LABELS = {"CV", "Max", "Mean", "Median", "SD"}

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def drop_leading_zeros(seg: pd.DataFrame, eps: float = 1e-6) -> pd.DataFrame:
    """Remove leading ~zero points (plotting only)."""
    vals = seg["Intensity"].to_numpy(dtype=float)
    mask = np.abs(vals) > eps
    if not mask.any():
        return seg.iloc[0:0].copy()
    first = int(np.argmax(mask))
    return seg.iloc[first:].copy()


def safe_name(text: str) -> str:
    """Filename-safe string."""
    bad = '[]:*?/\\'
    s = "".join(ch for ch in str(text) if ch not in bad)
    return s.replace(" ", "_")

def export_sample_map_template(samples: list[str], out_path: Path) -> None:
    """
    Write a user-editable mapping table: Sample -> (Amyloid, Dye)

    You fill in Amyloid and Dye, save, then re-run with --sample-map.
    """
    df = pd.DataFrame(
        {"Sample": sorted(set(str(s).strip() for s in samples)),
         "Amyloid": "",
         "Dye": ""}
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.suffix.lower() in [".xlsx", ".xls"]:
        with pd.ExcelWriter(out_path, engine="xlsxwriter") as w:
            df.to_excel(w, sheet_name="sample_map", index=False)
    else:
        df.to_csv(out_path, index=False)

    print(f"✅ Wrote sample map template: {out_path}")


def apply_sample_map(tidy: pd.DataFrame, map_path: Path) -> pd.DataFrame:
    """
    Merge a user-provided sample map into tidy and use it to set Amyloid/Dye.

    Behavior:
    - Keeps all tidy rows (left join on Sample)
    - If mapping provides a non-empty Amyloid/Dye, it overrides tidy's existing values
    - Warns if any samples are still missing labels
    """
    if map_path.suffix.lower() in [".xlsx", ".xls"]:
        m = pd.read_excel(map_path, sheet_name=0)
    else:
        m = pd.read_csv(map_path)

    required = {"Sample", "Amyloid", "Dye"}
    missing = required - set(m.columns)
    if missing:
        raise ValueError(f"Sample map missing columns: {sorted(missing)}")

    m = m[["Sample", "Amyloid", "Dye"]].copy()
    m["Sample"] = m["Sample"].astype(str).str.strip()

    # Treat blanks as NaN so they don't overwrite existing tidy values
    def _clean_col(col: pd.Series) -> pd.Series:
        col = col.astype(str).str.strip()
        col = col.replace({"": np.nan, "nan": np.nan, "None": np.nan})
        return col

    m["Amyloid"] = _clean_col(m["Amyloid"])
    m["Dye"] = _clean_col(m["Dye"])

    out = tidy.copy()
    out["Sample"] = out["Sample"].astype(str).str.strip()

    # Merge mapping columns as Amyloid_map/Dye_map
    out = out.merge(m.rename(columns={"Amyloid": "Amyloid_map", "Dye": "Dye_map"}),
                    on="Sample", how="left")

    # Override only where map has a value
    if "Amyloid" in out.columns:
        out["Amyloid"] = out["Amyloid_map"].combine_first(out["Amyloid"])
    else:
        out["Amyloid"] = out["Amyloid_map"]

    if "Dye" in out.columns:
        out["Dye"] = out["Dye_map"].combine_first(out["Dye"])
    else:
        out["Dye"] = out["Dye_map"]

    out = out.drop(columns=["Amyloid_map", "Dye_map"])

    # Warn about unmapped
    missing_mask = out["Dye"].isna() | (out["Dye"].astype(str).str.strip() == "") | \
                   out["Amyloid"].isna() | (out["Amyloid"].astype(str).str.strip() == "")
    if missing_mask.any():
        unmapped = out.loc[missing_mask, "Sample"].dropna().unique()
        print(f"⚠️ {len(unmapped)} samples still missing Amyloid and/or Dye.")
        print("   Examples:", list(unmapped[:10]))
        print("   (These samples may be skipped or grouped incorrectly in dye-based plots.)")

    return out


def make_plot_dirs(base_name: str, outdir: Path) -> Tuple[Path, Dict[str, Path]]:
    """Mirror spectral_preprocessing.py directory layout."""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    root = outdir / f"{ts}_{base_name}_sawtooth_plots"

    subdirs = {
        "uncollapsed_raw": root / "uncollapsed_raw",
        "uncollapsed_norm": root / "uncollapsed_norm",
        "collapsed_raw": root / "collapsed_raw",
        "collapsed_globalnorm": root / "collapsed_globalnorm_stacked",
        "collapsed_staggered": root / "collapsed_staggered_norm",
    }
    for p in subdirs.values():
        p.mkdir(parents=True, exist_ok=True)

    return root, subdirs


def _safe_token_case(s: str) -> str:
    """Title-case only if fully lowercase; otherwise keep original."""
    s = (s or "").strip()
    if not s:
        return s
    return s.title() if s.islower() else s


def parse_sample_fields(sample: str) -> Tuple[str, str, str]:
    """Parse Sample into (SampleClean, Amyloid, Dye) using '[Amyloid-Dye]' convention."""
    sample_clean = str(sample).strip().strip("[]")
    amyloid, dye = sample_clean, ""

    if "-" in sample_clean:
        left, right = sample_clean.split("-", 1)
        amyloid, dye = left.strip(), right.strip()

    return sample_clean, _safe_token_case(amyloid), _safe_token_case(dye)


def extract_channel_info(col_name: str) -> Optional[Tuple[str, int, float, int, str]]:
    """
    If col_name contains a channel token, return:
      (laser_prefix, channel_num, emission_start_nm, excitation_nm, measure (A, H, W))
    """
    m = CHANNEL_TOKEN_RE.search(col_name)
    if not m:
        return None

    prefix = m.group("prefix")
    num = int(m.group("num"))
    suffix = (m.group("suffix") or "")  # e.g. "-A", "-H", "-W"
    key = f"{prefix}{num}"
    if key not in CHANNEL_NM:
        return None

    emission_start = float(CHANNEL_NM[key])
    excitation = int(LASER_EXCITATION_NM[prefix])

    # normalize suffix to just "A/H/W" (or "" if none)
    measure = suffix[1:] if suffix.startswith("-") else suffix
    return prefix, num, emission_start, excitation, measure


# ---------------------------------------------------------------------
# NEW: wide-column renaming WITH excitation-row support
# ---------------------------------------------------------------------
def rename_wide_columns_with_excitation(
    cols: List[str],
    unique_names: bool = True,
) -> Tuple[List[str], List[float]]:
    """
    Rename channel tokens in wide headers and also return Excitation_nm values per column.

    Returns
    -------
    new_cols : list[str]
        Column names with channel tokens replaced by emission start (nm)
    excitation_vals : list[float]
        For each column in new_cols:
          - excitation nm if it was a channel column
          - NaN otherwise
    """
    new_cols: List[str] = []
    excitation_vals: List[float] = []

    for c in cols:
        ex_val = np.nan

        def _sub(match: re.Match) -> str:
            nonlocal ex_val
            prefix = match.group("prefix")
            num = int(match.group("num"))
            suffix = match.group("suffix") or ""
            key = f"{prefix}{num}"
            if key not in CHANNEL_NM:
                return match.group(0)

            start_nm = int(CHANNEL_NM[key])
            ex = LASER_EXCITATION_NM[prefix]
            ex_val = float(ex)

            # Either "498-A" or "488ex_498-A"
            return f"{ex}ex_{start_nm}{suffix}" if unique_names else f"{start_nm}{suffix}"

        new_c = CHANNEL_TOKEN_RE.sub(_sub, c)
        new_cols.append(new_c)
        excitation_vals.append(ex_val)

    return new_cols, excitation_vals


def write_cleaned_wide_csv(
    df: pd.DataFrame,
    out_csv: Path,
    unique_names: bool = False,
    include_excitation_row: bool = True,
) -> None:
    """
    Write a cleaned wide CSV:
      - rename channel tokens -> emission start (nm)
      - optionally prepend a first data row with excitation nm values
    """
    new_cols, excitation_vals = rename_wide_columns_with_excitation(
        list(df.columns), unique_names=unique_names
    )

    wide = df.copy()
    wide.columns = new_cols

    if include_excitation_row and len(new_cols) > 0:
        # Put a label in the first column for readability.
        first_col = new_cols[0]

        excitation_row = {col: val for col, val in zip(new_cols, excitation_vals)}
        excitation_row[first_col] = "Excitation_nm"

        wide = pd.concat([pd.DataFrame([excitation_row]), wide], ignore_index=True)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(out_csv, index=False)


# ---------------------------------------------------------------------
# Tidy building + normalization
# ---------------------------------------------------------------------
def wide_to_tidy(df: pd.DataFrame) -> pd.DataFrame:
    """Convert one wide Cytek Aurora dataframe into tidy rows for plotting."""
    if "Sample" in df.columns:
        sample_col = "Sample"
    elif "File Name" in df.columns:
        sample_col = "File Name"
    else:
        sample_col = df.columns[0]

    channel_cols = []
    for c in df.columns:
        info = extract_channel_info(c)
        if info is None:
            continue
        measure = info[-1]
        if measure != "A":     # <-- keep only Area
            continue
        channel_cols.append(c)

    if not channel_cols:
        return pd.DataFrame(
            columns=[
                "Sample", "SampleClean", "Amyloid", "Dye",
                "Excitation_nm", "Detector_CH", "Emission_nm_start",
                "Intensity", "OrderInBlock", "ChannelIndex"
            ]
        )

    records = []
    for _, row in df.iterrows():
        sample = str(row[sample_col]).strip()
        if sample in BAD_LABELS:
            continue

        sample_clean, amyloid, dye = parse_sample_fields(sample)

        for col in channel_cols:
            info = extract_channel_info(col)
            if info is None:
                continue
            _, detector_ch, emission_start, excitation, measure = info

            v = row.get(col, np.nan)
            try:
                intensity = float(str(v).replace(",", ""))
            except Exception:
                intensity = np.nan

            records.append(
                {
                    "Sample": sample,
                    "SampleClean": sample_clean,
                    "Amyloid": amyloid,
                    "Dye": dye,
                    "Excitation_nm": excitation,
                    "Detector_CH": int(detector_ch),
                    "Emission_nm_start": float(emission_start),
                    "Intensity": intensity,
                }
            )

    tidy = pd.DataFrame.from_records(records)
    if tidy.empty:
        return tidy

    tidy = tidy.sort_values(
        ["Sample", "Excitation_nm", "Emission_nm_start", "Detector_CH"],
        ascending=[True, True, True, True],
    ).reset_index(drop=True)

    tidy["OrderInBlock"] = tidy.groupby(["Sample", "Excitation_nm"]).cumcount() + 1
    tidy["ChannelIndex"] = tidy["OrderInBlock"]
    return tidy


def add_normalized_intensity(tidy: pd.DataFrame) -> pd.DataFrame:
    """
    Match spectral_preprocessing.py:
      normalize within each (Amyloid, Dye, Excitation_nm) group.
    """
    tidy = tidy.copy()

    def _norm(series: pd.Series) -> pd.Series:
        m = series.max(skipna=True)
        if not np.isfinite(m) or m <= 0:
            return pd.Series([0.0] * len(series), index=series.index)
        return series / m

    tidy["Intensity_norm"] = (
        tidy.groupby(["Amyloid", "Dye", "Excitation_nm"])["Intensity"]
            .transform(_norm)
    )
    return tidy


# ---------------------------------------------------------------------
# Plotting (same as previous enhanced version)
# ---------------------------------------------------------------------
def _get_excitation_order(tidy: pd.DataFrame) -> List[int]:
    present = set(int(x) for x in tidy["Excitation_nm"].dropna().unique())
    ordered = [ex for ex in PREFERRED_EXCITATION_ORDER if ex in present]
    extras = sorted(present - set(ordered))
    return ordered + extras


def plot_sawtooths_per_dye(tidy: pd.DataFrame, dye_name: str, outdir: Path, normalize: bool = False) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    sub = tidy[tidy["Dye"] == dye_name].copy()
    if sub.empty:
        return

    amyloids = sorted(sub["Amyloid"].dropna().unique())
    excitations_here = _get_excitation_order(sub)
    if not amyloids or not excitations_here:
        return

    n_rows, n_cols = len(amyloids), len(excitations_here)
    cmap = cm.get_cmap("viridis", len(excitations_here))
    ex_to_color = {ex: cmap(i) for i, ex in enumerate(excitations_here)}

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.2 * n_cols, 3.0 * n_rows),
        sharex=True, sharey=True
    )
    if n_rows == 1:
        axes = np.expand_dims(axes, axis=0)

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
            if normalize and np.nanmax(y) > 0:
                y = y / np.nanmax(y)

            ax.plot(x, y, marker="o", linewidth=2, color=ex_to_color[ex])

            if i == 0:
                ax.set_title(f"{ex} nm")
            if j == 0:
                ylabel = amy + ("\n(norm. max = 1)" if normalize else "")
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(axis="y", left=False, labelleft=False)

    title = f"Uncollapsed Sawtooths – {dye_name}" + (" (Normalized)" if normalize else "")
    fig.suptitle(title, fontsize=16, y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    outdir.mkdir(parents=True, exist_ok=True)
    fname = f"uncollapsed_sawtooths_{safe_name(dye_name)}{'_norm' if normalize else ''}.png"
    plt.savefig(outdir / fname, dpi=300)
    plt.close(fig)


def plot_collapsed_spectra_raw(tidy: pd.DataFrame, dye_name: str, outdir: Path) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    sub = tidy[tidy["Dye"] == dye_name].copy()
    if sub.empty:
        return

    excitations_here = _get_excitation_order(sub)
    amyloids = sorted(sub["Amyloid"].dropna().unique())
    if not amyloids or not excitations_here:
        return

    cmap = cm.get_cmap("viridis", len(excitations_here))
    outdir.mkdir(parents=True, exist_ok=True)

    for amy in amyloids:
        sub_amy = sub[sub["Amyloid"] == amy].copy()
        if sub_amy.empty:
            continue

        fig, ax = plt.subplots(figsize=(5, 4))
        plotted = False

        for i, ex in enumerate(excitations_here):
            seg = sub_amy[sub_amy["Excitation_nm"] == ex].copy()
            if seg.empty:
                continue
            seg = seg.sort_values("Emission_nm_start")
            seg = drop_leading_zeros(seg)
            if seg.empty:
                continue

            x = seg["Emission_nm_start"].to_numpy(dtype=float)
            y = seg["Intensity"].to_numpy(dtype=float)
            ax.plot(x, y, marker="o", linewidth=2, color=cmap(i), label=f"{ex} nm")
            plotted = True

        if not plotted:
            plt.close(fig)
            continue

        ax.set_xlabel("Emission (nm)")
        ax.set_ylabel("Intensity (A.U.)")
        ax.set_title(f"{amy} – {dye_name} (collapsed, raw)")
        ax.legend(title="Excitation nm", fontsize=8)
        ax.grid(False)
        plt.tight_layout()

        fname = f"collapsed_raw_{safe_name(dye_name)}_{safe_name(amy)}.png"
        plt.savefig(outdir / fname, dpi=300)
        plt.close(fig)


def plot_collapsed_spectra_global_norm(tidy: pd.DataFrame, dye_name: str, outdir: Path, offset_step: float = 0.5) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    sub = tidy[tidy["Dye"] == dye_name].copy()
    if sub.empty:
        return

    excitations_here = _get_excitation_order(sub)
    amyloids = sorted(sub["Amyloid"].dropna().unique())
    if not amyloids or not excitations_here:
        return

    cmap = cm.get_cmap("viridis", len(excitations_here))
    outdir.mkdir(parents=True, exist_ok=True)

    for amy in amyloids:
        sub_amy = sub[sub["Amyloid"] == amy].copy()
        if sub_amy.empty:
            continue

        fig, ax = plt.subplots(figsize=(5, 4))
        max_y_seen = 0.0
        plotted = False

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

            curve_max = np.nanmax(y_raw)
            if not np.isfinite(curve_max) or curve_max <= 0:
                continue

            y_norm = y_raw / curve_max
            y = y_norm + i * offset_step

            max_y_seen = max(max_y_seen, np.nanmax(y))
            ax.plot(x, y, marker="o", linewidth=2, color=cmap(i), label=f"{ex} nm")
            plotted = True

        if not plotted:
            plt.close(fig)
            continue

        ax.set_xlabel("Emission (nm)")
        ax.set_ylabel("Normalized (per excitation) + offset")
        ax.set_ylim(0, max_y_seen + 0.1)
        ax.set_title(f"{amy} – {dye_name} (collapsed, per-curve norm + stacked)")
        ax.legend(title="Excitation nm", fontsize=8)
        ax.grid(False)
        plt.tight_layout()

        fname = f"collapsed_globalnorm_{safe_name(dye_name)}_{safe_name(amy)}.png"
        plt.savefig(outdir / fname, dpi=300)
        plt.close(fig)


def plot_collapsed_spectra_staggered_norm(tidy: pd.DataFrame, dye_name: str, outdir: Path) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    sub = tidy[tidy["Dye"] == dye_name].copy()
    if sub.empty:
        return

    excitations_here = _get_excitation_order(sub)
    k = len(excitations_here)
    amyloids = sorted(sub["Amyloid"].dropna().unique())
    if not amyloids or k == 0:
        return

    cmap = cm.get_cmap("viridis", k)

    if k == 1:
        target_maxima = [1.0]
    else:
        step = 0.9 / (k - 1)
        target_maxima = [1.0 - i * step for i in range(k)]

    outdir.mkdir(parents=True, exist_ok=True)

    for amy in amyloids:
        sub_amy = sub[sub["Amyloid"] == amy].copy()
        if sub_amy.empty:
            continue

        fig, ax = plt.subplots(figsize=(5, 4))
        plotted = False

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

            curve_max = np.nanmax(y_raw)
            if not np.isfinite(curve_max) or curve_max <= 0:
                continue

            desired_max = target_maxima[i]
            y = (y_raw / curve_max) * desired_max

            ax.plot(x, y, marker="o", linewidth=2, color=cmap(i), label=f"{ex} nm (max={desired_max:.2f})")
            plotted = True

        if not plotted:
            plt.close(fig)
            continue

        ax.set_xlabel("Emission (nm)")
        ax.set_ylabel("Staggered normalized intensity (A.U.)")
        ax.set_ylim(bottom=0)
        ax.set_title(f"{amy} – {dye_name} (collapsed, staggered norm)")
        ax.legend(title="Excitation nm", fontsize=7)
        ax.grid(False)
        plt.tight_layout()

        fname = f"collapsed_staggered_{safe_name(dye_name)}_{safe_name(amy)}.png"
        plt.savefig(outdir / fname, dpi=300)
        plt.close(fig)


# ---------------------------------------------------------------------
# Excel exports
# ---------------------------------------------------------------------
def export_tidy_excel(tidy: pd.DataFrame, excel_path: Path) -> None:
    excel_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        tidy.to_excel(writer, sheet_name="tidy", index=False)


def export_sawtooths_by_dye_excel(tidy: pd.DataFrame, excel_path: Path) -> None:
    dyes = sorted(tidy["Dye"].dropna().unique())
    if not dyes:
        return

    excel_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        for dye in dyes:
            sub = tidy[tidy["Dye"] == dye].copy()
            if sub.empty:
                continue

            sub = sub.sort_values(["Amyloid", "Excitation_nm", "Emission_nm_start"])
            cols = [
                "Sample", "Amyloid", "Dye", "Excitation_nm",
                "Emission_nm_start", "Detector_CH", "OrderInBlock",
                "Intensity", "Intensity_norm",
            ]
            cols = [c for c in cols if c in sub.columns]
            sub = sub[cols]

            safe = "".join(ch for ch in str(dye) if ch not in '[]:*?/\\')
            sheet_name = (safe[:31] or "Sheet")
            sub.to_excel(writer, sheet_name=sheet_name, index=False)


# ---------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------
def collect_inputs(inputs: List[str], input_folder: Optional[str]) -> List[Path]:
    paths: List[Path] = [Path(p) for p in (inputs or [])]
    if input_folder:
        folder = Path(input_folder)
        if folder.is_dir():
            paths.extend(sorted(folder.glob("*.csv")))
        else:
            raise SystemExit(f"--input-folder not found or not a directory: {folder}")

    seen = set()
    out: List[Path] = []
    for p in paths:
        p = p.resolve()
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def process_and_combine(csv_paths: List[Path]) -> pd.DataFrame:
    tidies = []
    for csv_path in csv_paths:
        df = pd.read_csv(csv_path, low_memory=False)
        tidy = wide_to_tidy(df)
        if not tidy.empty:
            tidies.append(tidy)
    return pd.concat(tidies, ignore_index=True) if tidies else pd.DataFrame()


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Convert Cytek Aurora wide CSV(s) to cleaned wide CSVs + tidy Excel + sawtooth plots."
    )
    ap.add_argument("--inputs", nargs="*", default=[], help="One or more CSV files.")
    ap.add_argument("--input-folder", default=None, help="Folder containing CSVs (*.csv).")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--name", default=None, help="Base name for combined outputs (default: stem or 'combined').")

    ap.add_argument("--no-plots", action="store_true", help="Skip plot generation.")
    ap.add_argument("--unique-names", action="store_true", help="Embed excitation in wide renamed headers.")

    # NEW: manual labeling workflow
    ap.add_argument(
        "--export-sample-map",
        default=None,
        help="Write a Sample->(Amyloid,Dye) template (csv or xlsx) and exit."
    )
    ap.add_argument(
        "--sample-map",
        default=None,
        help="Path to a filled-in sample map (csv or xlsx) used to label Amyloid/Dye before plotting."
    )

    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    csv_paths = collect_inputs(args.inputs, args.input_folder)
    if not csv_paths:
        raise SystemExit("No input CSVs found. Provide --inputs and/or --input-folder.")

    # (Optional) your cleaned-wide export step stays here if you want it
    # for csv_path in csv_paths:
    #     df = pd.read_csv(csv_path, low_memory=False)
    #     write_cleaned_wide_csv(...)

    # Build combined tidy from inputs
    tidy = process_and_combine(csv_paths)
    if tidy.empty:
        raise SystemExit("No channel columns found (no UV/V/B/YG/R tokens detected).")

    # ===== NEW STEP A: Export editable sample map and exit =====
    if args.export_sample_map:
        map_path = Path(args.export_sample_map)
        export_sample_map_template(
            samples=tidy["Sample"].dropna().unique().tolist(),
            out_path=map_path
        )
        print("✋ Fill in Amyloid/Dye in the exported sample map, then re-run with --sample-map.")
        return

    # ===== NEW STEP B: Apply manual labels if provided =====
    if args.sample_map:
        tidy = apply_sample_map(tidy, Path(args.sample_map))

    # Now proceed as usual: normalize, export, plot
    tidy_norm = add_normalized_intensity(tidy)

    base_name = args.name or (csv_paths[0].stem if len(csv_paths) == 1 else "combined")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    tidy_norm_path = outdir / f"{ts}_{base_name}_tidy_norm.xlsx"
    export_tidy_excel(tidy_norm, tidy_norm_path)

    dye_excel_path = outdir / f"{ts}_{base_name}_sawtooths_by_dye_norm.xlsx"
    export_sawtooths_by_dye_excel(tidy_norm, dye_excel_path)

    if not args.no_plots:
        plot_root, plot_dirs = make_plot_dirs(base_name, outdir)
        print(f"📁 Writing plots to: {plot_root}")

        # IMPORTANT: dye list comes from tidy_norm AFTER mapping
        dyes = sorted(tidy_norm["Dye"].dropna().unique())
        dyes = [d for d in dyes if str(d).strip() != "" and str(d).lower() != "nan"]

        for dye in dyes:
            plot_sawtooths_per_dye(tidy_norm, dye, plot_dirs["uncollapsed_raw"], normalize=False)
            plot_sawtooths_per_dye(tidy_norm, dye, plot_dirs["uncollapsed_norm"], normalize=True)

            plot_collapsed_spectra_raw(tidy_norm, dye, plot_dirs["collapsed_raw"])
            plot_collapsed_spectra_global_norm(tidy_norm, dye, plot_dirs["collapsed_globalnorm"])
            plot_collapsed_spectra_staggered_norm(tidy_norm, dye, plot_dirs["collapsed_staggered"])

    print(f"✅ Wrote tidy Excel: {tidy_norm_path}")
    print(f"✅ Wrote by-dye Excel: {dye_excel_path}")
    print("Done.")


if __name__ == "__main__":
    main()
