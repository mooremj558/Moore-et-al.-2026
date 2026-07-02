file=$(awk -v  awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' settings.txt)

python /wanglab/mmoore/myapps/bin/3DNetMod_GPS_MMCP.py $file
