#!/bin/bash
#SBATCH --job-name=hicluster_generatematrix
#SBATCH --array=1-4%4
#SBATCH --time=12-0
#SBATCH --ntasks=1
#SBATCH -c 2
#SBATCH --mem=80G
##SBATCH --output=hicluster_impute100kb_%A-%a.out

#Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cLoops2

cLoops2_callLoops_array.sh
