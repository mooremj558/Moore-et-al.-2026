#!/bin/bash
#SBATCH --job-name=cLoops2_diffloop
#SBATCH --array=2-10%5
#SBATCH --time=12-0
#SBATCH --ntasks=1
#SBATCH -c 10
#SBATCH --mem=200G
#SBATCH --partition=gpu,general
##SBATCH --output=cLoops2_quantloop_%A-%a.out

#Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cLoops2

cLoops2_quantLoops_array.sh
