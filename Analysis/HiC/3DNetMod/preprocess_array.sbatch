#!/bin/bash
#SBATCH --job-name=3DNetMod_pre
#SBATCH --array=1-276,369-460%10
#SBATCH --time=12-0
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=gpu,general
##SBATCH --output=3DNetMod_pre_%A-%a.out

#Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cremins

preprocess_array.sh
