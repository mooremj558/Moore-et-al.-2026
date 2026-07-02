#!/bin/bash
#SBATCH --job-name=cLoops2_diffloop
#SBATCH --array=1-1%1
#SBATCH --time=12-0
#SBATCH --ntasks=1
#SBATCH -c 10
#SBATCH --mem=450G
#SBATCH --partition=gpu,general
##SBATCH --output=cLoops2_diffloop_%A-%a.out

#Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cLoops2

cLoops2_diffLoops_array.sh
