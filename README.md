# EMBER Dyes Amyloid Discrimination Analysis

## 🌟 Overview

This repository contains a complete Python pipeline for **processing fluorescence spectral data** from amyloid–binding dye experiments.  
The workflow transforms raw instrument CSV outputs into analyzed visualizations and statistical scores, enabling comparisons between amyloid strains and dyes — similar to the **EMBER (Yang et al., 2023, PNAS)** workflow.

## 🧰 Installation

### 1️⃣ Clone the repository

```bash
git clone https://github.com/<your-username>/Spectral_Data_Processing.git
cd Spectral_Data_Processing
```

### 2️⃣ Set up a virtual environment (recommended)

macOS/Linux
python3 -m venv venv
source venv/bin/activate

Windows
python -m venv venv
venv\Scripts\activate

### 3️⃣ Install dependencies
pip install -r requirements.txt

🚀 Step-by-Step Usage
Step 1 — Process Raw ClarioStar Spectral CSVs

Script: spectral_processing_clarioStar.py

Processes raw plate-reader CSVs from the ClarioStar and extracts emission–excitation spectral data.

python spectral_processing_clarioStar.py


Output:
Cleaned run-level CSV files (e.g., Spillover_LSM_A_processed.csv)

Step 2 — Preprocess and Normalize Spectra

Script: spectral_preprocessing.py

Cleans, reshapes, and normalizes the processed spectra into a tidy format suitable for statistical analysis.
Generates per-dye “sawtooth” plots for visual inspection.

python spectral_preprocessing.py


Output files:

combined_samples_tidy.xlsx

combined_samples_tidy_norm.xlsx

Per-dye Excel sheets with spectra (*_sawtooths_by_dye.xlsx)

PNG sawtooth plots in a _sawtooth_plots directory

Step 3 — Combine Multiple CSVs or Folders

Script: spectral_combine_csv.py

Merges data from multiple runs or folders into a single combined dataset.
Automatically detects sample rows formatted as [Amyloid-Dye].

Example 1: Use default paths (defined in script)
python spectral_combine_csv.py

Example 2: Specify folders or CSVs from terminal
python spectral_combine_csv.py /path/to/run1 /path/to/run2 /path/to/file3.csv


Each path should be separated by a space.
If a path contains spaces, wrap it in quotes.

Output:
combined_samples.csv containing all [amyloid-dye] rows with SourceFile and SourceFolder metadata.

Step 4 — Perform PCA, UMAP, and QDA Analysis

Script: spectral_analysis.py

Loads a tidy Excel file (from Step 2 or 3) and performs:

Principal Component Analysis (PCA)

UMAP dimensionality reduction

Quadratic Discriminant Analysis (QDA) for class separability

Supports multiple replicates per amyloid (e.g., asyn1-bf188, asyn2-bf188),
and automatically creates:

AmyloidBase → base amyloid group (e.g., “asyn”)

RunNumber → replicate index (e.g., “1”, “2”)

Example
python spectral_analysis.py combined_samples_tidy.xlsx


Output:

PCA and UMAP scatterplots (color-coded by amyloid)

QDA discrimination scores printed per dye, e.g.:

=== bf188 ===
PCA-QDA discrimination score: 0.872
  asyn: 0.910
  tau: 0.845

🧠 Notes and Best Practices

Sample names must follow the [Amyloid-Dye] pattern (e.g. [asyn-bf188]), with no underscores.

Multiple runs for the same amyloid should be labeled like asyn1-bf188, asyn2-bf188, etc.
The code automatically groups them by amyloid base for statistical analysis.

Always inspect generated sawtooth plots to ensure spectra are correctly aligned before analysis.

When running scripts from the terminal, separate input paths with spaces.
