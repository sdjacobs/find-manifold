#!/bin/bash

#SBATCH --job-name=manifold
#SBATCH --array=1000,2000,4000,8000
#SBATCH --output=bench/opt-%A_%a.txt
#SBATCH --error=bench/opt-%A_%a.err
#SBATCH --partition=bigmem
#SBATCH --nodes=1

module load python/2.7-2014q2

/usr/bin/time python main.py --output outtest $SLURM_ARRAY_TASK_ID 9 11
