"""
Script for running swat
To set it up ---> Go to setup.py
To submit as job to hpc ---> Go to submitHPC.py

Jose P. Teran - Nov 2025
"""


from setup import *
import os
import argparse
import multiprocessing
from datetime import datetime
import subprocess

SCRIPT_DIR = os.path.dirname(__file__)
LOG_DIR = os.path.join(SCRIPT_DIR, "logs")
os.makedirs(LOG_DIR, exist_ok=True)

def run_coswatv2(region):

    # per-region log folder
    date_str = datetime.now().strftime("%Y%m%d")
    region_log_dir = os.path.join(LOG_DIR,region,f"CoSWATv{version}")
    os.makedirs(region_log_dir, exist_ok=True)
    log_path = os.path.join(region_log_dir, f"{region}_{date_str}_out")

    with open(log_path, "w") as f:
        f.write(f"=== Region {region} started at {datetime.now()} ===\n")
        f.flush()

        cmds = [["python3", "write_swatplus_files.py", "--r", region]]

        if copy_new_files:
            cmds.append(["python3", "copyncheck_weather.py", "--r", region])

        cmds.append(["python3", "run_sims.py", "--r", region])

        for cmd in cmds:
            f.write(f"\n$ {' '.join(cmd)}\n")
            f.flush()

            # all prints + errors from the script go into the same log
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=False)

        f.write(f"\nRegion {region} finished at {datetime.now()}\n")



if all_regions:
    print('\n Running CoSWAT+ v2 for ALL regions ... \n\n')

    with open(all_regions_list_file,'r') as f:
        regions_raw = f.readlines()
        regions = []
        for region in regions_raw:
            region_ok = region.split('\n')[0]
            regions.append(region_ok)


else:
    print('\n Running CoSWAT+ v2 for seleced regions ... \n\n')
    regions = selected_regions


processes = nr_processes


if __name__ == "__main__":

    with multiprocessing.Pool(int(processes)) as pool:
        pool.map(run_coswatv2, regions)