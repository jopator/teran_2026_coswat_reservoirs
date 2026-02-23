# Integrating reservoirs and lakes in the CoSWAT Global Hydrological Model

This study enhances the CoSWAT Global Water Model by improving the representation of reservoirs and lakes at large scale. We implemented operational schemes for storage and release, and developed a topology-based approach to estimate irrigation demand for irrigation reservoirs using global datasets. The model was restructured and evaluated across multiple regions and compared the results with other global models.

------------------------------------------------------------------------

The datasets required to run the scripts can be obtained from this [Zenodo Repository](https://doi.org/10.5281/zenodo.18733431) (Teran,2026).

This study also used the [CoSWAT-Framework](https://github.com/jopator/CoSWAT-Framework). To fully reproduce the model setup from scratch, this framework is
required.

Moreover, a modified version of the Soil and Water Assessment Tool (SWAT+) was used.The source code can be accessed in [this repository](https://github.com/jopator/swatplus).

------------------------------------------------------------------------

# Reproduction Workflow

The following documents provide step-by-step instructions to reproduce
the results of this study starting from scratch:

-   [Setting up regions with the CoSWAT-Framework](docs/Setting_up_CoswatFramework.md)
-   [Running CoSWAT simulations](docs/running_sims.md)
-   [Post-processing CoSWAT outputs](docs/post_process_outputs.md)

The following document provides instructions to reproduce the analyses
presented in the manuscript using the post-processed outputs provided in
[Teran(2026)](https://doi.org/10.5281/zenodo.18733431) . This can be done without performing the full model setup and simulation
steps above.

-   [Analysis of results](docs/analysis_of_results.md)

------------------------------------------------------------------------

# General Instructions to Run Scripts

Please ensure that the folder structure is preserved.\
Files downloaded from [Teran(2026)](https://doi.org/10.5281/zenodo.18733431) must maintain their
relative structure and be merged appropriately with this repository.

To run a script from the command line:

``` bash
cd Scripts/result_analysis/reservoirs_coswat/validation
python3 00_run_parallel_val.py
```

Some scripts require arguments. These are described at the beginning of
each script.\
Typically, a CoSWAT region name and model version must be specified:

``` bash
cd Scripts/result_analysis/reservoirs_coswat/validation
python3 01_reservoir_storage_validation.py --r <region> --v <version>
```

------------------------------------------------------------------------

# Regions Used in This Study

-   africa-nile
-   africa-orange
-   america-bravo
-   america-colorado
-   america-mississippi
-   america-parana
-   asia-mekong
-   europe-central
-   europe-west

------------------------------------------------------------------------

# Model Versions Used in This Study

**1.5.0**\
Model setup at 0.01° resolution including all developments presented in
this study (lakes, reservoirs, and irrigation topology).

**1.1.0**\
Model setup at 0.01° resolution without the developments introduced in
this study.

------------------------------------------------------------------------

# Installing the Python Environment

A requirements.txt file is provided to install the Python packages used in a virtual environment.

To create the environment:

``` bash
python3 -m venv <virtual_environment_name>
source <virtual_environment_name>/bin/activate
pip install -r requirements.txt
```

Ensure that Python (3.10 or higher) is properly installed before executing the
commands above.

# Author
- Jose P. Teran
- Email: jose.pablo.teran.orsini@vub.be