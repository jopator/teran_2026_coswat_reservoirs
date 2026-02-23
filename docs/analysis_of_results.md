# Validation and Analysis

This section describes the scripts used to validate CoSWAT simulations and to perform performance analyses. These scripts operate on the processed outputs generated in the previous step.

All scripts are located in:

`Scripts/result_analysis/reservoirs_coswat/`

Validation requires either (i) locally processed simulation outputs or (ii) the processed outputs provided in [Teran(2026)](https://doi.org/10.5281/zenodo.18733431).

---

# Pre-processing Reservoir Reference Data (Optional)

Location:

`Scripts/res_obs_preprocessing/pre-process-reservoirs.py`

This script prepares the reference reservoir datasets used for validation. It harmonizes and aggregates data from:

* [GRS](https://doi.org/10.5281/zenodo.7855477)
* [ResOpsUS](https://doi.org/10.5281/zenodo.6612040)
* [Yassin (2018)](https://doi.org/10.5281/zenodo.1492043)

It filters out locations based on data availability (at least 1 year of continuous data and 5 years with at least one data point), and it also does a consistency evaluation, comparing long-term storage with reported capacity from the GranD dataset, filtering out reservoirs with large discrepancies.

These data are already processed and provided in [Teran (2026)](https://doi.org/10.5281/zenodo.18733431). Therefore, running this script is **not required** to reproduce the manuscript results. It is provided to allow full replication of the data preparation workflow.

To execute:

```bash
cd Scripts/res_obs_preprocessing
python3 pre-process-reservoirs.py
```

---

# Model performance

Location:

`Scripts/result_analysis/reservoirs_coswat/validation`

These scripts evaluate model performance using processed simulation outputs.

---

## 00_run_parallel_val.py

Utility script to execute validation routines in parallel across multiple regions.

Example:

```bash
python3 00_run_parallel_val.py
```

---

## 01_reservoir_storage_validation.py

Evaluates reservoir and lake performance.

Calculates:

* KGE
* PBIAS
* Correlation (R)

For:

* Storage
* Inflow
* Outflow

Uses processed outputs and the prepared reference datasets.

Arguments are `--r <region> --v <version>`

Example:

```bash
python3 01_reservoir_storage_validation.py --r africa-nile --v 1.5.0
```

Outputs:

* CSV files containing performance metrics per water body
* Time series comparison plots

---

## 02_streamflow_validation.py

Evaluates channel streamflow performance.

Calculates:

* KGE
* PBIAS
* Correlation (R)

Uses processed outputs that already include matched observation data.

Arguments are `--r <region> --v <version>`

Example:

```bash
python3 02_streamflow_validation.py --r africa-nile --v 1.5.0
```

Outputs:

* CSV files with performance metrics per station
* Time series plots

---

# Reservoir Performance Analysis

### 03_reservoir_performance_analysis.py

Aggregates and summarizes validation results for reservoirs and lakes [Figure 5, Figure B-1, Figure B-2].

Generates:

* Performance distribution plots
* Summary statistics
* Regional comparison figures

```bash
python3 03_reservoir_performance_analysis.py
```

---

### 04_reservoir_time_series_plot.py

Generates time series plots for selected reservoirs. It includes storage, inflow and outflow [Figure 6].

```bash
python3 04_reservoir_time_series_plot.py
```

---

# Streamflow Performance Analysis

### 05_streamflow_performance_processing.py

Filters out channel stations influenced by reservoirs prior to performance aggregation.

Example:

```bash
python3 05_streamflow_performance_processing.py
```

---

### 06_streamflow_performance_analysis.py

Generates summary performance plots and distributions of metrics [Figure 8]. It generates a plot of streamflow climatologies for selected stations [Figure 9].

```bash
python3 06_streamflow_performance_analysis.py
```

---

# ISIMIP Model Comparison

Location:

`Scripts/result_analysis/reservoirs_coswat/validation/isimipModels`

This directory contains scripts to evaluate outputs from five ISIMIP global models:

* CWATM
* H08
* LPJmL5-7-10-fire
* MIROC-INTEG-LAND
* WaterGAP2-2e

The model outputs can be obtained at the [ISIMIP3A Global Water Sector repository](https://doi.org/10.48364/ISIMIP.398165.9). Use the "reservoirstor" variable, and download the files of each model to `.../validation/isimipModels/modelOutputs`.

---

## processOutputs.py

Samples ISIMIP model outputs at the same locations used for CoSWAT validation and prepares a pickle file with processed data per CoSWAT region at `.../validation/isimipModels/ProcessedData`.

```bash
cd isimipModels/evaluation
python3 processOutputs.py
```

---

## validationISIMIP.py

Computes performance metrics (KGE, PBIAS, R) for the ISIMIP models at the same validation locations and creates .csv files with the results at `.../validation/isimipModels/stats`.

```bash
python3 validationISIMIP.py
```

---

## summary_analysis.py
Generates a summary plot of the model performance evaluation [Figure 7], and generates a log file with a detailed report.

```bash
python3 summary_analysis.py
```

---

# Additional Analyses

Location:

`Scripts/result_analysis/reservoirs_coswat/others`

These scripts perform supplementary analyses in the study.

---

## resolve_efficiency.py

Evaluates the efficiency of resolving lakes and reservoirs into the river network.

```bash
python3 resolve_efficiency.py
```

---

## irrigation_demand_nasser.py

Estimates irrigation demand applied downstream of Lake Nasser in the Nile River Basin.

```bash
python3 irrigation_demand_nasser.py
```

---

# Helper Script

### read_swat.py

Contains helper classes and functions for reading SWAT+ output files.
Imported internally by the validation and processing scripts.

This script should not be executed independently.