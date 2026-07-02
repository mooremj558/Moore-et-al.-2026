#!/bin/bash

#SBATCH --time=240:00:00
#SBATCH -J juicer_cleanup
#SBATCH --partition=general

module load juicer

#/opt/apps/juicer/1.6.2/scripts/cleanup.sh
/scratch/devtools/yliang/softwares/juicer/scripts/cleanup.sh
