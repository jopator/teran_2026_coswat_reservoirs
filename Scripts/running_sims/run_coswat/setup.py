"""
Run settings setup to run coswatv2

Jose P. Teran - Nov 2025
"""

from pathlib import Path
import os

BASE_DIR = Path(__file__).resolve().parents[3]  # root of repo
os.chdir(BASE_DIR)

print(BASE_DIR)
# General Settings
version             = '1.5.0'
exe_path            = 'Scripts/running_sims/exes/swatplus-61.0.2.11-coswatv1.5-ifx-lin_x86_64' #'/data/brussel/vo/000/bvo00033/vsc10883/PHD_main/Projects/CoSWAT/Scripts/running_sims/exes/swatplus-61.0.2.11-coswatv1.5-ifx-lin_x86_64' #'/vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00033/CoSWAT-Framework/data-preparation/resources/swatplus-61.0.2-CoSWATv0.5.0-rel-ifx-lin'
sim_regions_path    = f'CoSWAT-Framework/model-setup/CoSWATv{version}'

time_sim    = {
    'yrc_start'     : '1965',
    'yrc_end'       : '2015',
}

print_prt   = {
    'nyskip'        :'5',
    'csvout'        :'n',
    'basin_wb'      : {'daily':'n','monthly':'n','yearly':'y','avann':'y'},
    'channel_sd'    : {'daily':'n','monthly':'y','yearly':'y','avann':'y'},
    'channel'       : {'daily':'n','monthly':'y','yearly':'y','avann':'y'},
    'reservoir'     : {'daily':'n','monthly':'y','yearly':'y','avann':'y'}
}

codes_bsn   = {
    'pet'       : '1',
    'rte_cha'   : '0'
}

nr_processes = 9

# Simulation regions
all_regions             = False
all_regions_list_file   = 'Scripts/running_sims/regions.txt' #'/data/brussel/vo/000/bvo00033/vsc10883/PHD_main/Projects/CoSWAT/regions.txt'

selected_regions        = [
    'africa-nile',
    'africa-orange',
    'america-bravo',
    'america-colorado',
    'america-mississippi',
    'america-parana',
    'asia-mekong',
    'europe-central',
    'europe-west'
]

# Weather settings
copy_new_files          = False
source_dir              = '/data/brussel/vo/000/bvo00033/data_repo/weather/coswatv2/gswp3-w5e5' # Must specify if need to copy (if made with CoSWAT-Framework this is not needed)
path_to_files           = '{region}' #{region}/counterclim/20CRv3


# HPC Job settings
file_name   = f'Scripts/running_sims/run_coswat/runCoSWATv{version}sh'
partition   = 'zen5_himem'
nodes       = 1
ntasks      = nr_processes
job_name    = "runCoSWAT"
time        = "0-8:00:00"
output      = f"coswat_{version}"
