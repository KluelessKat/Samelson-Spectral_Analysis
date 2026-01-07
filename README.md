# EMBER Dyes Amyloid Discrimination Analysis

## 🌟 Overview

This repository contains a complete Python pipeline for **processing fluorescence spectral data** from amyloid–binding dye experiments.  
The workflow transforms raw instrument CSV outputs into analyzed visualizations and statistical scores, enabling comparisons between amyloid strains and dyes — similar to the **EMBER (Yang et al., 2023, PNAS)** workflow.

## 🧰 Installation

### 1️⃣ Clone the repository

```bash
git clone https://github.com/KluelessKat/Samelson-Spectral_Analysis.git
cd Samelson-Spectral_Analysis
```

### 2️⃣ Set up a virtual environment (recommended)

macOS/Linux
```bash
python3 -m venv venv
source venv/bin/activate
```

Windows
```bash
python -m venv venv
venv\Scripts\activate
```

### 3️⃣ Install dependencies
```bash
pip install -r requirements.txt
```

## 🚀 Step-by-Step Usage for SONY ID7000 Data
### Step 1 — Process Raw SONY ID7000 Spectral CSVs

Script: spectral_combine_csv.py
* The unmixed data from ID7000 has many rows that provide no information about our samples, so this script combines just the useful sample data from multiple unmixing .csv files and outputs one singular .csv file of the combined raw data.
* Merges data from multiple runs or folders into a single combined dataset.
* Automatically detects sample rows formatted as [Amyloid-Dye] (not case-sensitive).

Example 1: Use default paths (defined in script)
```bash
python spectral_combine_csv.py
```

Example 2: Specify folders or CSVs from terminal
```bash
python spectral_combine_csv.py /path/to/run1 /path/to/run2 /path/to/file3.csv
```
* Each path should be separated by a space.
* If a path contains spaces, wrap it in quotes.

Output:
* Cleaned run-level CSV file (e.g., combined_samples.csv) containing all [amyloid-dye] rows with SourceFile and SourceFolder metadata.

### Step 2 — Preprocess, Normalize Spectra, and Plot Sawtooth Curves

Script: spectral_preprocessing.py
* Cleans, reshapes, and normalizes the processed spectra into a tidy format suitable for statistical analysis.
* Generates per-dye “sawtooth” plots for visual inspection.
* Takes in processed csv output from spectral_combine_csv.py as the argument

```bash
python spectral_preprocessing.py /path/to/*_combined_samples.csv
```

Output files:
* Two per-dye Excel sheets with spectra (*_sawtooths_by_dye.xlsx)
  * *_combined_samples_tidy_norm.xlsx
  * *_combined_samples_tidy_norm.xlsx 
* PNG sawtooth plots in a *_sawtooth_plots directory

### Step 3 — Perform PCA, UMAP, and QDA Analysis

Script: spectral_analysis.py

* Loads a tidy Excel file (from Step 2) and performs:

* Principal Component Analysis (PCA)

* UMAP dimensionality reduction

* Quadratic Discriminant Analysis (QDA) for class separability

* Supports multiple replicates per amyloid (e.g., asyn1-bf188, asyn2-bf188),
and automatically creates:

* AmyloidBase → base amyloid group (e.g., “asyn”)

* RunNumber → replicate index (e.g., “1”, “2”)

Example
```bash
python spectral_analysis.py combined_samples_tidy.xlsx
```

Output:

PCA and UMAP scatterplots (color-coded by amyloid)
  * within directory named *_analysis_plots

QDA discrimination scores printed per dye, e.g.:

=== bf188 ===
PCA-QDA discrimination score: 0.872
  asyn: 0.910
  tau: 0.845

## 🧠 Notes and Best Practices

Sample names must follow the [Amyloid-Dye] pattern (e.g. [asyn-bf188]).

Multiple runs for the same amyloid should be labeled like asyn1-bf188, asyn2-bf188, etc.
The code automatically groups them by amyloid base for statistical analysis.

Always inspect generated sawtooth plots to ensure spectra are correctly aligned before analysis.

When running scripts from the terminal, separate input paths with spaces.
---
## 🚀 Step-by-Step Usage for Cytek Aurora Data
### Step 1 — Process Raw Cytek Aurora Spectral CSVs Exported from FloJo11

Script: spectral_ca_preprocessing.py

Script converts exported tables from FloJo11 to:
- A **tidy** (long-format) xlsx file for analysis and plotting
- A **normalized** tidy xlsx file (per excitation curve)
- **Sawtooth / collapsed spectrum plots** (per dye and amyloid)
- Excel exports (one combined tidy sheet, where there's one workbook with one sheet per dye)

The script is designed to be used **before** any downstream embedding / PCA / UMAP / classification steps.

### 1) What the script expects as input

### Supported input format
You provide one or more CSVs exported from your Cytek Aurora workflow where:

- One column contains the **sample name** (usually `Sample` or `File Name`).
- The spectral measurement columns include **Aurora channel tokens** such as:
  - `UV1-A`, `UV2-A`, …
  - `V1-A`, `B3-A`, `YG7-A`, `R8-A`, etc.

Only columns that match the pattern below are treated as channels:

- Prefix: `UV`, `V`, `B`, `YG`, `R`
- Channel number: `1–16` (varies by prefix)
- Optional suffix: typically `-A`, `-H`, `-W`

**Important:** the script only keeps **Area** measurements (`-A`). Height/Width (`-H`, `-W`) are ignored.

### Sample naming convention (automatic labeling)
By default, the script tries to infer:
- `Amyloid`
- `Dye`

from the `Sample` string using the convention:

```
[Amyloid-Dye]
```

Examples:
- `[Beads_2-Bf188]` → Amyloid = `Beads_2`, Dye = `Bf188`
- `[AD-THT]` → Amyloid = `AD`, Dye = `THT`

If your sample names do **not** follow this convention, use the **sample map** workflow (recommended).

---

### 2) Quickstart

### A) Run on a folder of CSVs (plots + Excel outputs)

```bash
python spectral_ca_preprocessing.py \
  --input-folder ./csvs \
  --outdir ./out
```

### B) Run on explicit files

```bash
python spectral_ca_preprocessing.py \
  --inputs ./csvs/run1.csv ./csvs/run2.csv \
  --outdir ./out \
  --name my_run
```
OR (See Step 5) Create Sample Map for Amyloid and Dye Labeling per Sample (Recommend)
```bash
python spectral_ca_preprocessing.py \
  --inputs ./csvs/run1.csv ./csvs/run2.csv \
  --outdir ./out \
  --name my_run \
  --export-sample-map ./edited_csv.csv
```

Then:
--> Manually add into edited_csv.csv (or whatever file name) the amyloid and dye
--> Save changes
--> Run same script again but with argument --sample-map as follows:

```bash
python spectral_ca_preprocessing.py \
  --inputs ./csvs/run1.csv ./csvs/run2.csv \
  --outdir ./out \
  --name my_run \
  --sample-map ./edited_csv.csv
```


### C) Skip plot generation (Excel only)

```bash
python spectral_ca_preprocessing.py \
  --input-folder ./csvs \
  --outdir ./out \
  --no-plots
```

---

### 3) Command-line arguments (full reference)

### Input selection
- `--inputs <file1.csv file2.csv ...>`  
  Provide one or more CSV files.

- `--input-folder <path>`  
  Provide a folder containing CSVs (`*.csv`). All CSVs in the folder are processed.

> You can use **either** `--inputs` **or** `--input-folder` (or both). Duplicate paths are automatically removed.

### Output configuration
- `--outdir <path>` (**required**)  
  Directory where all outputs will be written.

- `--name <string>`  
  Base name used in output filenames.
  - If you provide **one input file** and no `--name`, the script uses the file’s stem.
  - If you provide **multiple input files** and no `--name`, the script uses `combined`.

### Plotting
- `--no-plots`  
  If set, the script will not generate PNG sawtooth plots (takes a while to render).

### Header renaming (optional / advanced)
- `--unique-names`  
  Used for “cleaned wide CSV export” mode: embed excitation into renamed channel headers (e.g., `488ex_508-A`) so that wavelengths shared across laser families don’t collide.
  
> Note: using unique-names is currently set as True by default. This can be changed under the main() function.

### Sample labeling workflow (recommended)
- `--export-sample-map <path.(csv|xlsx)>`  
  Creates a **template** mapping file with columns: `Sample`, `Amyloid`, `Dye`, then exits. This allows you to manually add the amyloid and dye for each sample without worrying about reformatting the sample name to adhere to the "Amyloid-Dye" convention required by the script.

- `--sample-map <path.(csv|xlsx)>`  
  Applies your filled-in mapping file to label samples before plotting and exporting.

---

### 4) Outputs

All outputs are written under `--outdir`.

### A) Combined tidy Excel (`*_tidy_norm.xlsx`)
A single Excel file containing one sheet named `tidy` with one row per:

> (Sample × Excitation × Channel)

Columns include (at minimum):

- `Sample` (original label)
- `SampleClean` (brackets removed)
- `Amyloid`
- `Dye`
- `Excitation_nm`
- `Detector_CH` (channel number within laser family)
- `Emission_nm_start` (see note below)
- `Intensity`
- `Intensity_norm` (see normalization section)
- `OrderInBlock`, `ChannelIndex` (ordering helpers)

**Emission wavelength note:** the current script uses the **channel center wavelength** lookup table for `Emission_nm_start`. If you want start wavelengths instead, see Section 9.

### B) By-dye Excel (`*_sawtooths_by_dye_norm.xlsx`)
One workbook containing one sheet per dye, each listing all rows for that dye (sorted by amyloid, excitation, emission).

### C) Plot directory (`*_sawtooth_plots/`)
Unless you pass `--no-plots`, the script creates a timestamped plot folder with subfolders:

- `uncollapsed_raw/`  
- `uncollapsed_norm/`  
- `collapsed_raw/`  
- `collapsed_globalnorm_stacked/`  
- `collapsed_staggered_norm/`  

Each dye gets:
- “Uncollapsed” sawtooths: grid with **rows = amyloid**, **cols = excitation**
- “Collapsed” spectra: one figure per (amyloid, dye), with multiple excitations overlaid

---

### 5) Sample map workflow (highly recommended)

If your sample naming is inconsistent, or if you want to control grouping/labels exactly, use a sample map.

### Step 1 — export a template
```bash
python spectral_ca_preprocessing.py \
  --input-folder ./csvs \
  --outdir ./out \
  --export-sample-map ./out/sample_map.csv
```

This creates a file with columns:
- `Sample` (auto-filled, unique values found in your CSVs)
- `Amyloid` (blank)
- `Dye` (blank)

### Step 2 — fill in labels
Open `sample_map.xlsx` and fill `Amyloid` and `Dye` for each `Sample`.

Rules:
- Leave cells blank if you do **not** want to override the auto-parsed label (when sample name is already in `Amyloid-Dye` format).
- You may assign multiple Samples to the same Amyloid/Dye to group replicates.
- Don’t rename the headers: the script requires `Sample`, `Amyloid`, `Dye`.

### Step 3 — run using the mapping file
```bash
python spectral_ca_preprocessing.py \
  --input-folder ./csvs \
  --outdir ./out \
  --sample-map ./out/sample_map.csv
```

The merge behavior is:

- Left-join on `Sample` (no tidy rows are dropped)
- If `Amyloid`/`Dye` is present in the map, it overrides the parsed value
- The script prints a warning if any samples still have missing `Amyloid`/`Dye` after mapping

---

### 6) Optional: exporting “cleaned wide CSVs” (advanced)

The script includes helper functions to export a “cleaned wide CSV” where channel headers are renamed to emission wavelength and an optional first row lists excitation wavelengths.

However, this portion is currently commented out in `main()`.

If you want this behavior, you can uncomment and implement the following block in `main()`:

```python
# for csv_path in csv_paths:
#     df = pd.read_csv(csv_path, low_memory=False)
#     write_cleaned_wide_csv(...)
```

### Why `--unique-names` matters for wide CSV export
If you rename headers to just the emission wavelength (e.g. `508-A`), you may get **duplicate column names** across laser families (UV/V/B/YG/R can share wavelength bins). Embedding excitation avoids collisions:

- Without unique names: `508-A` may appear multiple times → columns collide
- With unique names: `488ex_508-A`, `405ex_508-A`, etc. → unique

---
### 8) Run PCA/UMAP/QDA Analysis
script: spectral_analysis.py

```bash
python spectral_ca_preprocessing.py ./out/combined_samples_tidy_norm.xlsx
```

### 8) Troubleshooting / common issues

### “No channel columns found (no UV/V/B/YG/R tokens detected).”
This means your CSV headers did not match the expected token pattern (e.g., `UV1-A`).

Check:
- Are you exporting the right table from Aurora/FlowJo?
- Do the channel columns include prefixes like `UV`, `V`, `B`, `YG`, `R`?
- Do the columns include the `-A` suffix? (The script filters to Area only.)

### Dyes are missing / everything groups under one dye
If your `Sample` strings don’t contain a `-` separator inside brackets (e.g., `[Amyloid-Dye]`), the script can’t infer `Dye`.

Fix:
- Use the sample map workflow (Section 7).

### Plots appear to “drop leading points”
The plotting code intentionally removes leading near-zero intensities (for readability). This does **not** modify the exported Excel values; it only affects plots.

### Some wavelengths seem to “repeat” or “disappear” in wide-format outputs
This typically happens when multiple channels map to the same numeric wavelength and you rename wide headers without embedding excitation.

Fix:
- Use excitation-embedded names (`--unique-names`) for wide header renaming, or avoid wide renaming and rely on tidy format instead.




