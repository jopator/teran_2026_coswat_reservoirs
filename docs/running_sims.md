# Running CoSWAT Simulations

This section describes how to execute CoSWAT simulations after the model setup has been completed (see Model Setup section).

Only two actions are required:

1. Configure the simulation in `setup.py`
2. Execute `coswatv2.py` (or submit the provided `.sh` as slurm jobs if on an HPC system)

All other scripts in this directory are internal components of the workflow and do not need to be executed manually.

---

### Directory Structure
    running_sims
    ├── exes
    │ └── swatplus-61.0.2.11-coswatv1.5-ifx-lin_x86_64
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

## 1. Configure the Simulation

Before running simulations, you need to edit `run_coswat/setup.py`. The following parameters must be defined:

- **Model version** (`1.5.0` or `1.1.0`) at line 15.
- **Parallel processes** (Ideally 9) at line 38.


---

## 2. Running Simulations

To run simulations (locally) from the command line:

```bash
cd running_sims/run_coswat
python3 coswatv2.py
```

This script:
- Reads configuration from setup.py
- Iterates over selected regions
- Updates necessary SWAT+ input files
- Calls the SWAT+ executable on each region

No other Python scripts in this directory need to be executed manually.

In the case of HPC systems where a slurm job can be submitted, use:
- `runCoSWATv1.1.0.sh`
- `runCoSWATv1.5.0.sh`

This scripts internally execute `coswatv2.py`.

### SWAT+ Executable

The executable used in this study is located in:

`running_sims/exes/swatplus-61.0.2.11-coswatv1.5-ifx-lin_x86_64`

It was compiled for Linux systems using cmake and the intel fortran compiler and includes the modifications introduced in this study.

To compile the code yourself (if needed), please follow the guide provided in the [repository](https://github.com/jopator/swatplus).

---

### Internal Scripts

The following scripts are called automatically by coswatv2.py and should not be executed individually:

* write_swatplus_files.py
* copyncheck_weather.py
* jopato_fx.py
* run_sims.py
* create_bash.py

### Expected Output

If execution is successful, simulation outputs will be written to:

`model-setup/CoSWATv<version>/<region>/Scenarios/<scenario>/TxtInOut/`

These outputs are used in the post-processing and analysis steps described in the subsequent documentation.
