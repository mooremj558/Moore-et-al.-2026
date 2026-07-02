#!/bin/bash
#SBATCH --job-name=3DNetMod_GPS
#SBATCH --array=415-552%10
#SBATCH --time=12-0
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --partition=gpu,general
##SBATCH --output=3DNetMod_GPS_%A-%a.out

#Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cremins

HSVM_array.sh
