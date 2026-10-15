# Laundry Wastewater and River-System Water Analysis

This project analyzes water samples made available by participating laundries and collected in connection with a river and locations along its course. Organic and inorganic substance tests are used to characterize the samples and support future wastewater-treatment planning.

The tracked dataset is organized by laundry and sample condition. It contains untreated wastewater (`ST`), treated water (`T`), and sludge (`L`) samples. The analysis compares chemical detection patterns across these sample types and laundries. It is an analytical study; the repository does not document a designed treatment process.

## Data

`muestras.csv` contains laboratory results for 97 compounds:

- 29 PFAS (compound IDs 1–29)
- 68 restricted substances (compound IDs 30–97)

The current file includes 13 sample columns from five laundries: 7 untreated wastewater samples, 5 treated-water samples, and 1 sludge sample. Sample IDs identify laundry and sample type; geographic coordinates for individual collection sites are not included in the tracked dataset.

Results below the reporting threshold are marked with `<` and treated as censored observations in the scripts.

## Analysis workflow

1. Read compound and sample metadata from `muestras.csv`.
2. Summarize detections by compound, laundry, and sample type.
3. Compare chemical profiles and detection consistency across laundries.
4. Select 24 compounds detected in at least three laundries for focused analysis.
5. Export summary tables and generate detection, profile, histogram, and boxplot figures.

`eda.py` uses the threshold divided by two for censored values in its concentration matrix. `imputacion.py` applies the US EPA substitution approach implemented in the script, using the threshold divided by √2 for the selected untreated-water dataset. These substitutions are estimates, not measured concentrations. The EDA script notes substantial censoring and a small untreated-sample count, so imputed results should be treated as exploratory rather than definitive.

## Repository contents

- `muestras.csv` — source laboratory results and sample labels
- `eda.py` — exploratory analysis, detection summaries, filtered datasets, and figures
- `imputacion.py` — censored-value imputation and concentration-distribution figures
- `figuras/` — generated and saved visualizations
- `outputs/` — filtered datasets and exported summaries

## Reproduce the exploratory analysis

Run these commands from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install numpy pandas matplotlib seaborn
python eda.py
```

The EDA script reads `muestras.csv` and writes plots to `figuras/` and tables to `outputs/`.

To run the imputation workflow, first generate `outputs/muestras_24_st.csv` with `eda.py`. The current `imputacion.py` expects a file named `muestras_24_st.csv` in its working directory, while the EDA script writes that file under `outputs/`; place the generated file at the expected path before running:

```bash
python imputacion.py
```

The imputation workflow saves its dataset and additional figures under `outputs/` and `figuras/`.
