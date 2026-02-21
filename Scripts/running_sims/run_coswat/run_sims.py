"""
Run SWAT+ in parallel

Jose P. Teran - Nov 2025
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
import setup
import os
from jopato_fx import *
import argparse
from datetime import datetime
import coswatv2 as settings


# set default working directory to the script location

if __name__ == '__main__':

    # argument and regions
    parser = argparse.ArgumentParser(description="scripts for running coswatv2")
    parser.add_argument("--r", help="region", type=str, default='')
    args = parser.parse_args()

    if len(args.r) > 0: region = args.r
    else: quit()

    date_str = datetime.now().strftime("%Y%m%d")

    SCRIPT_DIR = os.path.dirname(__file__)
    LOG_DIR = os.path.join(SCRIPT_DIR, "logs")
    region_log_dir = os.path.join(LOG_DIR,region,f"CoSWATv{setup.version}")
    log_path = os.path.join(region_log_dir, f"{region}_{date_str}_out")

    # variables
    workingDir  = setup.sim_regions_path                                                                            
    exePath     = setup.exe_path                                                                                 
    txtInOut = f"{region}/Scenarios/Default/TxtInOut"

    run_dir = f"{workingDir}/{txtInOut}"

    runModel(run_dir,exePath,log_path)
