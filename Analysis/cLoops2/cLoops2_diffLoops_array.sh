file=$(awk -v  awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' celllines.txt)

cLoops2 callDiffLoops -tloop /scratch/mmoore/Epitherapy3D/analysis/HiC/4_loop/cloops_chr/merged_loops/${file}_DMSO_loops_CTCF.txt -cloop /scratch/mmoore/Epitherapy3D/analysis/HiC/4_loop/cloops_chr/merged_loops/${file}_DP_loops_CTCF.txt -td /scratch/mmoore/Epitherapy3D/analysis/HiC/4_loop/cloops/${file}_DMSO/ -cd /scratch/mmoore/Epitherapy3D/analysis/HiC/4_loop/cloops/${file}_DP/ -o diffloops/new_unsampled/${file}_DMSOvsDP -p 10 -mcut 10000000 -vmin="-1" -vmax=3
