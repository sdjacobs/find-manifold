#!/bin/bash

#SBATCH --job-name=manifold
#SBATCH --output=singlecore.txt
#SBATCH --error=singlecore.err
#SBATCH --partition=sandyb
#SBATCH --nodes=1

module load python/2.7-2014q2

python findManifold.py -m 2000 9
