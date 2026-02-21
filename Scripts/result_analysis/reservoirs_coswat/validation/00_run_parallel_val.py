"""
Author: Jose P. Teran
Date: 2025-11-28

Description:
Run in parallel the validation scripts for selected regions
"""

import os
import multiprocessing
from datetime import datetime
import subprocess


# Setup
regions = ['africa-nile',
           'africa-orange',
           'america-bravo',
           'america-colorado',
           'america-mississippi',
           'america-parana',
           'asia-mekong',
           'europe-central',
           'europe-west'
           ]

version   = "1.5.0"
processes = 9

object_validation = "reservoir" # or "channel"

SCRIPT_DIR = os.path.dirname(__file__)
LOG_DIR = os.path.join(SCRIPT_DIR, "logs")
os.makedirs(LOG_DIR, exist_ok=True)

if object_validation == "reservoir":
    python_script = "01_reservoir_storage_validation.py"
else:
    python_script = "02_streamflow_validation.py"

def run_validation(region):
    # per-region log folder
    date_str = datetime.now().strftime("%Y%m%d")
    region_log_dir = os.path.join(LOG_DIR,f"CoSWATv{version}",f"{region}")
    os.makedirs(region_log_dir, exist_ok=True)
    log_path = os.path.join(region_log_dir, f"{region}_{date_str}_out")

    with open(log_path, "w") as f:

        cmds = [
            ["python3", os.path.join(SCRIPT_DIR, python_script),
             "--r", region, "--v", version]
            ]

        for cmd in cmds:
            f.write(f"\n$ {' '.join(cmd)}\n")
            f.flush()

            # all prints + errors from the script go into the same log
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=False)

        f.write(f"\nRegion {region} finished at {datetime.now()}\n")

if __name__ == "__main__":

    with multiprocessing.Pool(int(processes)) as pool:
        pool.map(run_validation, regions)