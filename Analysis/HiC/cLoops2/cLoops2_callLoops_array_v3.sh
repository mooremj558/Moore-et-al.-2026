file=$(awk -v  awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' samples2.txt)

cLoops2 callLoops -d ${file} -o ${file} -eps 5000,7500,10000 -minPts 5,20,50 -p 15 -w -j -hic -mcut 10000000 -max_cut
