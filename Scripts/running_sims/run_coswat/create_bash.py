"""
Create bash file

Jose P. Teran - Nov 2025
"""

import os
from setup import *

script = f"""#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --job-name={job_name}
#SBATCH --time={time}
#SBATCH --output={output}

module load Python/3.11.3-GCCcore-12.3.0
source /vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00033/vsc10883/PHD_main/jopato_env_3.11/bin/activate

echo "CoSWAT+ HPC Simulation Framework"
chmod +x ./coswatv2.py
python3 coswatv2.py
"""


with open(file_name, "w") as f:
    f.write(script)

os.chmod(file_name, 0o755)