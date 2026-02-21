#!/bin/bash
#SBATCH --partition=zen5_himem
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --job-name=processOut
#SBATCH --time=0-2:00:00
#SBATCH --output=process_adjusted_coswat

module load Python/3.11.3-GCCcore-12.3.0
source /vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00033/vsc10883/PHD_main/jopato_env_3.11/bin/activate

echo "Post Processing CoSWAT+ Outputs - Reservoir and Channels"
chmod +x ./run_process.py
python3 run_process.py
