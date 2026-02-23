# Setting up regions with the CoSWAT-Framework

This section describes the full model setup from raw inputs.
If you only wish to reproduce the manuscript results using processed outputs, see: [Analysis of results](analysis_of_results.md).


This tool requires docker (https://www.docker.com/). Please make sure it is installed in your system. This study was done with Docker 29.2.1 on Ubuntu 24.04.

Full simulation setup is computationally intensive and may require HPC resources. Processed outputs are provided to allow reproduction of manuscript results without rerunning simulations. Model files and simulation outputs require at least 8 TB of storage.

# Contents
- [Preparation](#preparation)
- [Initial settings](#initial-settings)
- [Prepare forcings](#prepare-forcings)
- [Model Setup](#model-setup)

## Preparation
- Download the files from the [repository](https://github.com/jopator/CoSWAT-Framework). In this repository, you can also find a general guide on its usage.

- Save those files, respecting the structure in the "CoSWAT-Framework" of this repository and including the data provided in [Teran (2026)](https://doi.org/10.5281/zenodo.18733431).

You should at least have the following folder structure:

    ├── CoSWAT-Framework
    │   ├── buildDocker
    │   ├── data-preparation
    │   ├── Dockerfile
    │   ├── docker-resources
    │   ├── main-scripts
    │   ├── model-data
    │   ├── README.md
    │   └── runDocker
    ├── Scripts
    │   ├── res_obs_preprocessing
    │   ├── result_analysis
    │   └── running_sims

**model-data** should have the contents downloaded from [Teran (2026)](https://doi.org/10.5281/zenodo.18733431) for the 9 regions.

- Build the docker image:
```bash
cd CoSWAT-Framework
chmod +x buildDocker
./buildDocker
```

## Initial Settings
The **datavariables.py** file, in the **main-scripts** directory is used to set up basic settings on the CoSWAT-Framework. For detailed information, please visit the official repository. Below, several settings associated with this study are explained.


#### Version Tag
In this variable you can provide a version tag to the setup. For this study, the tags are 1.5.0 for the version with reservoirs and lakes, and 1.1.0 without them.
```python
version  # line 11
```

#### Reservoir and lake filtering parameters

- grand_lake_thres: Area threshold to filter out small unregulated lakes.
- grand_min_thres : Area threshold to filter out small regulated reservoirs.
- grand_dor_thres : Threshold to filter out lakes and reservoirs with low degree of regulation.

```python
grand_lake_thres  # line 83
grand_min_thres   # line 84
grand_dor_thres   # line 85
```

#### Irrigation topology parameters

- main_threshold: Threshold downstream sub-basins in the irrigation topology.

- tributary_order: Tributary order level for streams considered on irrigation topology for a reservoir.

- fao_irrg_area_pctg_thres: Threshold for irrigated fraction in Global Irrigated Area Map to consider an agricultural HRU as irrigated.


```python
main_threshold   # line 93
tributary_order  # line 94
fao_irrg_area_pctg_thres  # line 95
```

#### Atmospheric forcing preparation
This study used the GSWP3-W5E5 dataset (Lange et al.,2022). For this, some specifications need to be done.

To allow the CoSWAT-Framework to download and process the atmospheric forcings in order to be used in coswat, these variables need to be set as follows:


```python
prepare_weather     = True # Line 106
redo_weather        = True # Line 107
weather_redownload  = True # Line 108

# ...

available_models    = ['gswp3-w5e5', ... # Line 121 

# ...


scenariosData = {                   # Line 131
    'observed'  : ['gswp3-w5e5',], 

    #...
}
```

## Prepare Forcings
The repository does not provide the atmospheric forcings directly in the SWAT+, so you need to use the CoSWAT-Framework to download the gridded files from the [ISIMIP repository](https://data.isimip.org/).

With the specified settings, atmospheric forcing data (GSWP3-W5E5) are automatically downloaded from the ISIMIP repository via the CoSWAT-Framework.

To do this, first, run a container with the created docker image:

```bash
chmod +x runDocker
./runDocker
```

Once in the container, run the necessary script:

```bash
cd data-preparation

python3 prepare-weather.py africa-nile africa-orange \
 america-bravo america-colorado america-mississippi america-parana \
 asia-mekong europe-central europe-west

```

This will create a **weather** folder with the processed data for each region in **model-data**.

## Model Setup
After the creation of weather files, you can follow with the model setup:
```bash
cd ../main-scripts

python3 set-up-model.py africa-nile africa-orange \
 america-bravo america-colorado america-mississippi america-parana \
 asia-mekong europe-central europe-west \
 --v <version>

```

The version should be either **1.5.0** or **1.1.0**.

If successful, the following folder should be created: `model-setup/CoSWATv1.5.0/<region>`.

For example, the america-bravo region, provided in the Zenodo Repository:

    CoSWATv1.5.0/america-bravo
    ├── america-bravo.qgs
    ├── america-bravo.qgs~
    ├── america-bravo.sqlite
    ├── DatabaseBackups
    ├── dem_data.xml
    ├── dem.qml
    ├── Scenarios
    ├── swatplus_datasets.sqlite
    └── Watershed
