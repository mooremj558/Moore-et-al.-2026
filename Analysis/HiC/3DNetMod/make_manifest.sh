#!/usr/bin/sh

for chr in {1..22} X
do
for plat in 16
do
for value in 4
do
for overlap in 400
do
for region_size in 600
do
for size in 15
do
for file in `ls *10k_VC_SQRT_interactions.txt`
do 
xbase=$(basename $file)
sample_temp=${xbase/_10k_VC_SQRT_interactions.txt/}
sample=${sample_temp/_/}
cat general_manifest.txt | sed "s/dummy_sample/${sample}/g" | sed "s/dummy_bed/bins_genomewide_10kb.bed/g" | sed "s/dummy_counts/${sample_temp}_10k_VC_SQRT_interactions.txt/g" | sed "s/dummy_chr/chr${chr}/g" | sed "s/dummy_plat/${plat}/g" | sed "s/dummy_overlap/${overlap}/g" | sed "s/dummy_region_size/${region_size}/g" | sed "s/dummy_size/${size}/g" > settings/settings_${sample}_chr${chr}_F_wonform_0.85_plat${plat}_value${value}_merge_${size}.txt
done
done
done
done
done
done
done
