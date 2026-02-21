# Description
This repository contains the scripts used for the study: Integrating reservoirs and lakes in the CoSWAT global hydrological model.\
This study also used the CoSWAT-Framework, which can be accessed at: https://github.com/jopator/CoSWAT-Framework
And the Soil and Water Assessment Tool (SWAT+), which can be accessed at: https://github.com/swat-model/swatplus

This study enhances the CoSWAT Global Water Model by improving the representation of reservoirs and lakes at large scale. We implemented operational schemes for storage and release, and developed a topology-based approach to estimate irrigation demand for irrigation reservoirs using global datasets. The model was restructured and evaluated across multiple regions and compared the results with other global models.


# General instructions to run scripts
- Please make sure to define references to other files correctly. Most scripts use the os.chdir() function to establish a working directory - All path definitions must be relative to this.
- Some scripts require arguments to be run individually. This is mentioned at the beginning of each script.
- A yaml file is provided to create a python environment suitable to run the scripts in case it is needed.

# Data availability.

To be able run simulations and execute processing/analysis scripts, the model files for this study are needed. They can be obtained in the following repository: https://doi.org/10.5281/zenodo.18633376 (Teran, 2026). Here, you can also find already processed outputs if you only wish to execute the analysis scripts without needing to run the model.

- The model files are provided in the "CoSWAT-Framework/model-setup" directory. They are available for CoSWAT version 1.5.0 (With lakes and reservoirs) and version 1.1.0(without lakes and reservoirs).

- Processed outputs are provided in the coswat_outputs.zip - These are processed outputs for selected waterbodies as part of the model evaluation section.

- Observations for reservoir storage, inflow and outflow are provided in the observations.zip. Observations for streamflow are provided in the "CoSWAT-Framework/model-data"

# PART 1: Running CoSWAT Simulations

The `running_sims` folder contains the tools required to configure and execute CoSWAT simulations across multiple regions.


## Folder Structure

    └── running_sims
        ├── exes
        │   └── swatplus-61.0.2.11-coswatv1.5-ifx-lin_x86_64
        ├── regions.txt
        └── run_coswat
            ├── copyncheck_weather.py
            ├── coswatv2.py
            ├── create_bash.py
            ├── jopato_fx.py
            ├── runCoSWATv1.1.0sh
            ├── runCoSWATv1.5.0sh
            ├── run_sims.py
            ├── setup.py
            └── write_swatplus_files.py
---

## Description

### `exes/`

Contains the compiled SWAT+ executable used in this study (swatplus-61.0.2.11-coswatv1.5-ifx-lin_x86_64`). It is a binary to be used in a Linux system.\
The source code can be obtained at the official SWAT+ repository provided above. The source code already includes additions made in this study, and can be compiled for other operative systems.

---

### `regions.txt`

Plain text file listing all CoSWAT regions to be simulated.  
---

## `run_coswat/`

Core scripts to configure and execute simulations.

#### `setup.py`

Defines simulation settings, including:
- CoSWAT version
- Regions to simulate
- General execution options

This file should be configured before launching simulations.

---

#### `coswatv2.py`

Main execution script.

It:
- Reads the configuration from `setup.py`
- Iterates over selected regions
- Calls the SWAT+ executable
- Manages the regional simulation workflow

This is the script that must be executed to run simulations.

---

#### `write_swatplus_files.py`

Helper script used to update basic SWAT+ configuration/input files within each regional folder prior to execution.

---

#### `copyncheck_weather.py`

Utility script to copy required weather input files from a reference source and check if they are properly structured.

---

#### `jopato_fx.py`

Collection of helper functions used by the execution scripts.

---

#### `run_sims.py`

Auxiliary script to trigger simulation workflows.

---

#### `runCoSWATv1.1.0sh`  
#### `runCoSWATv1.5.0sh`

Batch scripts for submission on HPC systems.\
These are used to launch CoSWAT simulations as scheduled jobs. They run coswatv2.py.


# PART 2: Processing outputs and performing analysis

This can be performed using the scripts and datasets associated with the results_analysis folder.


## Folder Structure

    Scripts/result_analysis/
    └── reservoirs_coswat
        ├── data_processing
        │   ├── process_out_coswat_reservoirs.py
        │   ├── process_out_coswat_streamflow.py
        │   ├── read_swat.py
        │   ├── run_process.py
        │   └── run_process.sh
        ├── others
        │   ├── irrigation_demand._nasser.ipynb
        │   ├── nasser_irrigated_hrus.gpkg
        │   ├── read_swat.py
        │   └── resolve_efficiency.ipynb
        └── validation
            ├── 00_run_parallel_val.py
            ├── 01_reservoir_storage_validation.py
            ├── 02_streamflow_validation.py
            ├── 03_reservoir_performance_analysis.py
            ├── 04_reservoir_time_series_plot.py
            ├── 05_streamflow_performance_processing.py
            ├── 06_streamflow_performance_analysis.py
            ├── read_swat.py
            ├── shapefiles
---

## `data_processing/`

Scripts to process outputs from CoSWAT simulations: they produce pickle files with time series per model object (i.e., reservoir or channel).

#### `process_out_coswat_reservoirs.py`
Processes out reservoir storage, inflow and outflow values. \
It takes --r <region> and --v <version> as argument

#### `process_out_coswat_streamflow.py`
Processes out flow values for selected channels and observations from the Global Runoff Data Centre (GRDC: https://grdc.bafg.de/data/data_portal/). \
This requires the CoSWAT-Framework folder provided in the associated repository.\
It takes --r <region> and --v <version> as argument

#### `run_process.py` & `run_process.sh`
This script can be used to run both processing scripts mentioned above in parallel for multiple regions.
The .sh file is a batch script for submission on HPC systems. It will run run_process.py.\

---
## `validation/`
Scripts to process outputs from CoSWAT simulations: they produce pickle files with time series per model object (i.e., reservoir or channel).

### `00_run_parallel_val.py`
This script can be use to run the validation scripts in parallel for multiple regions. They are described below.\

#### `01_reservoir_storage_validation.py`
This script calculates performance indices (KGE, PBIAS, R) and produces plots of the simulated and observed time series.\
This is done for reservoir/lake storage, inflow and outflow, using the processed outputs.\
This requires the observations folder with the dataset provided in the associated repository, which is an aggregation of the GRS(https://zenodo.org/records/7855477), ResOpsUS(https://zenodo.org/records/5893641) and Yassin 2018(https://zenodo.org/records/1492043) datasets.
It saves the performance indices per water body on a .csv file per region, as well as the plots.\
It takes --r <region> and --v <version> as argument

#### `02_reservoir_storage_validation.py`
This script calculates performance indices (KGE, PBIAS, R) and produces plots of the simulated and observed time series.\
This is done for streamflow, using the processed outputs.\
The processed outputs already include observations.\
It saves the performance indices per water body on a .csv file per region, as well as the plots.\
It takes --r <region> and --v <version> as argument

#### `03_reservoir_performance_analysis.py` & `04_reservoir_time_series_plot.py`
These scripts perform an analysis of the results of validation for reservoirs/lakes.\
The `03_reservoir_performance_analysis.py` generates a report and a number of plots associated with general performance for storage, inflow and outflow.
The `04_reservoir_time_series_plot.py` generates a plot of storage, inflow and outflow time series for selected reservoirs.

#### `05_streamflow_performance_processing.py` & `06_streamflow_performance_analysis.py`
These scripts perform an analysis of the results of validation for streamflow.
The `05_streamflow_performance_processing.py` filters out channels that are downstream of reservoirs.
The `06_streamflow_performance_analysis.py` generates report and a number of plots associated with general performance, as well as streamflow climatology for selected stations.

#### `read_swat.py`
This script has helper objects and functions to more easily read SWAT+ (hence CoSWAT) outputs. It is imported in many of the scripts mentioned above.

#### `shapefiles/`
This folder contains some vector files used for plotting. This includes river basins from HydroBASINS (https://www.hydrosheds.org/products/hydrobasins), major rivers in the world (https://datacatalog.worldbank.org/search/dataset/0042032/major-rivers-of-the-world) and lakes from HydroLAKES (https://www.hydrosheds.org/products/hydrolakes).