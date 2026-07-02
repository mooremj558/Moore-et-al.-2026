file=$(awk -v  awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' conditions.txt)

cLoops2 quant -d /scratch/mmoore/Epitherapy3D/analysis/HiC/4_loop/cloops/${file}/ -loops diffloops/3_merged_diffloop_quant/merged_diffloops.txt -o ${file} -p 10 -mcut 10000000

