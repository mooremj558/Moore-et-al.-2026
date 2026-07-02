file=$(awk -v  awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' settings_HSVM.txt)

python /wanglab/mmoore/myapps/bin/3DNetMod_HSVM_v2.py $file ${file/plat16/plat8}
