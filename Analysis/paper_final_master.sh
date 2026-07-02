###SVL analysis
cd /scratch/mmoore/Epitherapy3D/hic/paper_final/0_SVL/mega
mkdir ./VC_matrix
bash /scratch/yliang/Epitherapy_3D/OLD/scripts/6_1.dump_matrix.sh --maxjobs 10 --hicfolder /scratch/mmoore/Epitherapy3D/hic/paper_final/0_SVL/mega/hic --outputfolder /scratch/mmoore/Epitherapy3D/hic/paper_final/0_SVL/mega/VC_matrix --normalization VC_SQRT --oe observed --bin_size 50000 --bin_size_n 50k
rm VC_matrix/*chrY*
python /scratch/yliang/Epitherapy_3D/scripts/hic/0_QC/3_SVL/0_3.SVL_ratio.py --input /scratch/mmoore/Epitherapy3D/hic/paper_final/0_SVL/mega/VC_matrix --output 0_3.ALL

###SVL hicexplorer analysis
find ../hic -maxdepth 1 -name "*hic" | while read file; do xbase=$(basename $file); echo "hicConvertFormat -m $file  --inputFormat hic --outputFormat cool -o ${xbase/hic/cool} --resolutions 50000" >> hic2cool_commands.txt; done;
parallel_GNU -j 10 < hic2cool_commands.txt 2> hic2cool_commands.err &> hic2cool_commands.log
cd ../cool_50kb_VC
find ../cool_50kb -maxdepth 1 -name "*cool" | while read file; do xbase=$(basename $file); echo "hicConvertFormat -m $file  --inputFormat cool --outputFormat cool -o ${xbase/cool/VC.cool} --correction_name VC_SQRT" >> cool2cool_commands.txt ; done;
parallel_GNU -j 10 < cool2cool_commands.txt 2> cool2cool_commands.err &> cool2cool_commands.log
hicPlotSVL -m cool_50kb_VC/B36_DMSO_50000.VC.cool cool_50kb_VC/B36_DP_50000.VC.cool cool_50kb_VC/B49_DMSO_50000.VC.cool cool_50kb_VC/B49_DP_50000.VC.cool cool_50kb_VC/B66_DMSO_50000.VC.cool cool_50kb_VC/B66_DP_50000.VC.cool cool_50kb_VC/NHA_DMSO_50000.VC.cool cool_50kb_VC/NHA_DP_50000.VC.cool cool_50kb_VC/qNHA_DMSO_50000.VC.cool cool_50kb_VC/qNHA_DP_50000.VC.cool --distance 2000000 --threads 4 --plotFileName SVL_all.png --outFileName SVL_all_pval.txt --outFileNameData SVL_all_data.txt
hicPlotSVL -m cool_50kb_VC/B36_DMSO_50000.VC.cool cool_50kb_VC/B36_DP_50000.VC.cool cool_50kb_VC/B49_DMSO_50000.VC.cool cool_50kb_VC/B49_DP_50000.VC.cool cool_50kb_VC/B66_DMSO_50000.VC.cool cool_50kb_VC/B66_DP_50000.VC.cool cool_50kb_VC/NHA_DMSO_50000.VC.cool cool_50kb_VC/NHA_DP_50000.VC.cool cool_50kb_VC/qNHA_DMSO_50000.VC.cool cool_50kb_VC/qNHA_DP_50000.VC.cool --distance 50000 --threads 4 --plotFileName SVL_all_50.png --outFileName SVL_all_50_pval.txt --outFileNameData SVL_all_50_data.txt


###Compartment calling
cd /scratch/mmoore/Epitherapy3D/hic/paper_final/1_compartment/mega
for file in `ls *hic`; do xbase=$(basename $file ".hic"); echo "bash /bar/yliang/softwares/compartment_call_Bing_pipeline_2015/1_hic2mat.sh --sample $file --basename $xbase --matrix_path ./VC_matrix --bin_size 50000 --bin_size_n 50k --normalization VC_SQRT --REF_DIR /bar/yliang/genomes --juicer_path /bar/yliang/softwares/juicer/scripts/common --pipe_path /bar/yliang/softwares/compartment_call_Bing_pipeline_2015" >> hic_to_matrix_commands.txt; done;
mkdir ./VC_matrix
parallel_GNU -j 10 < hic_to_matrix_commands.txt 
bash /bar/yliang/softwares/compartment_call_Bing_pipeline_2015/HiC_pipe_compartment_splitarm_final.sh -j 10 -h /scratch/mmoore/Epitherapy3D/hic/paper_final/1_compartment/mega -m /scratch/mmoore/Epitherapy3D/hic/paper_final/1_compartment/mega/VC_matrix/ -b 50000 -bn 50k -n VC_SQRT -dm no

bash ../compartment_flip_atac.sh -b final.B36_DMSO_1.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B36_DMSO_BRep1_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B36_DMSO_2.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B36_DMSO_BRep2_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B36_DP_1.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B36_DP_BRep1_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B36_DP_2.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B36_DP_BRep2_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B49_DMSO_1.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B49_DMSO_BRep1_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B49_DMSO_2.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B49_DMSO_BRep2_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B49_DP_1.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B49_DP_BRep1_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B49_DP_2.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B49_DP_BRep2_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B66_DMSO_1.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B66_DMSO_BRep1_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B66_DMSO_2.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B66_DMSO_BRep2_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B66_DP_1.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B66_DP_BRep1_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B66_DP_2.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B66_DP_BRep2_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B66_DP_3.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B66_DP_BRep1_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw

bash ../compartment_flip_atac.sh -b final.B36_DMSO.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B36_DMSO_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B36_DP.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B36_DP_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B49_DMSO.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B49_DMSO_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B49_DP.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B49_DP_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B66_DMSO.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B66_DMSO_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.B66_DP.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_B66_DP_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.NHA_DMSO.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_NHA_DMSO_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.NHA_DP.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_NHA_DP_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.qNHA_DMSO.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_qNHA_DMSO_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw
bash ../compartment_flip_atac.sh -b final.qNHA_DP.bing2015.eigen.VC_SQRT.50k.bedgraph -a Trimmed_ATAC_qNHA_DP_R1.fastq_nochrM_nodup_qffilter_uniqpeak_treat_pileup_sorted.bw


###Compartment overlap with RNA/ATAC
for file in `ls ./atac_norm/` ; do echo $file ; xbase=$(basename $file) ; sample=${xbase/.bw/} ; bigWigAverageOverBed atac_norm/${file} 50kb_bins.bed ${sample}_comp_atac_avg.tab -bedOut=${sample}_comp_atac_avg.bed ; done

for file in `ls *DEGs.bed` ; do xbase=$(basename $file) ; tail -n +2 $file | bedtools intersect -a ../50kb_bins.bed -b stdin -wo > ${xbase/.bed/_compartment.bed} ; done

awk '{ $2=sprintf("%.0f",$2); $3=sprintf("%.0f",$3); print }' OFS="\t" /scratch/devtools/mmoore/genomes/private/hg38_CpG.bed | bedtools intersect -a ../atac_rna/50kb_bins.bed -b stdin -c > 50kb_bins_CpGnum.bed


###TAD calling
###Bing pipeline
for file in `ls *hic`; do xbase=$(basename $file ".hic"); echo "bash /bar/yliang/softwares/compartment_call_Bing_pipeline_2015/1_hic2mat.sh --sample $file --basename $xbase --matrix_path ./VC_matrix --bin_size 10000 --bin_size_n 10k --normalization VC_SQRT --REF_DIR /bar/yliang/genomes --juicer_path /bar/yliang/softwares/juicer/scripts/common --pipe_path /bar/yliang/softwares/compartment_call_Bing_pipeline_2015" >> hic_to_matrix_commands.txt; done;
parallel_GNU -j 10 < hic_to_matrix_commands.txt
bash /scratch/yliang/Epitherapy_3D/scripts/2_1.HiC_pipe_domain.sh --maxjobs 7 --hicpath /scratch/mmoore/Epitherapy3D/hic/paper_final/2_TAD/mega --dumpmatrix no --matrixpath /scratch/mmoore/Epitherapy3D/hic/paper_final/2_TAD/mega/VC_matrix --binsize 10000 --binsizen 10k --window_size_di 200 --window_size_ins 50 --sampletype our --normalization VC_SQRT

###3DNetMod
find ../hic -maxdepth 1 -name "*hic" | while read file; do xbase=$(basename $file); echo "hicConvertFormat -m $file  --inputFormat hic --outputFormat cool -o ${xbase/hic/cool} --resolutions 10000" >> hic2cool_commands.txt; done;
parallel_GNU -j 10 < hic2cool_commands.txt 2> hic2cool_commands.err &> hic2cool_commands.log
cd ../cool_10kb_VC
find ../cool_10kb -maxdepth 1 -name "*cool" | while read file; do xbase=$(basename $file); echo "hicConvertFormat -m $file  --inputFormat cool --outputFormat cool -o ${xbase/cool/VC.cool} --correction_name VC_SQRT" >> cool2cool_commands.txt ; done;
parallel_GNU -j 10 < cool2cool_commands.txt 2> cool2cool_commands.err &> cool2cool_commands.log
cd ../interactions
find ../cool_10kb_VC -maxdepth 1 -name "*cool" | while read file; do xbase=$(basename $file); echo "cooler dump -o ${xbase/10000.VC.cool/10k_VC_SQRT_interactions.txt} $file" >> coolerdump_commands.txt ; done;
parallel_GNU -j 10 < coolerdump_commands.txt 2> coolerdump_commands.err &> coolerdump_commands.log

bash make_manifest.sh 
bash make_manifest_reg300.sh 

cat settings/* > settings.txt
sbatch preprocess_array.sbatch
sbatch GPS_MMCP_array.sbatch
##HSVM at 7,70000; 10,100000; 15,150000 settings ; need to run NHA and qNHA in separate directories
sbatch HSVM_array.sbatch

for size in 70000 100000 150000 ; do for cell in B36 B49 B66 NHA qNHA ; do for chr in {1..22} X ; do for treat in DMSO DP ; do echo "cat ../output*/HSVM/*/merged/*chr${chr}${cell}${treat}*${size}*adjust3_b.txt >> "${cell}_${treat}_${size}.boundary ; done ; done ; done ; done > boundary_files_commands.txt
. boundary_files_commands.txt

for size in 70000 100000 150000 ; do for cell in B36 B49 B66 NHA qNHA ; do for chr in {1..22} X ; do for treat in DMSO DP ; do echo "cat ../output*/HSVM/*/merged/*chr${chr}${cell}${treat}*${size}*adjust3.txt >> "${cell}_${treat}_${size}.TAD ; done ; done ; done ; done > TAD_files_commands.txt
. TAD_files_commands.txt

#0_input folder with TAD files and DI files
cd 1_DI_quantile_normalization
#Quantile norm R notebook chunk
cd 2_boundaries_with_DI/
###Make merged raw TADs file
for file in ../0_input/*.TAD; do
  name=$(basename $file)
  echo $name 
  ln $file $name.bed
  done
## merge tad boundaries. 
#rm combined_tads.raw.sorted.txt
for file in $(ls *.TAD.bed); do xbase=$(basename $file)
  awk -v OFS="\t" -v name=${xbase/.TAD.bed/} '{ if(NR>1) print $0,name}' $file
    done |sort --parallel=4 -k1,1 -k2,2n -k3,3n - >> combined_tads.raw.sorted.txt
## Find replicated TAD boundaries 
for file in *.bed; do
file2=$(basename $file)
name=${file2/.TAD.bed/}
echo $name
awk -v OFS="\t" -v name=$name '{print $1,$2,name"\n"$1,$3,name}' $file |sort -k1,1 -k2,2n > ./$name.boundary 
done
sort -k1,1 -k2,2n -m ./*.boundary > ./boundary.all.txt
cat boundary.all.txt | uniq > boundary.all.uniq.txt
#Merge boundaries R notebook chunk
# overlap TAD boundary with DI score.
for file in $(ls ../1_DI_quantile_normalization/*.bdg); do xbase=$(basename $file) ; sample=${xbase/.normalized.DI.bdg/} ;
  echo $sample
  intersectBed \
    -a <(awk -v OFS="\t" '{if ($2-50000 <0){print $1,0,$2+50000,$1":"$2 } else {print $1,$2-50000,$2+50000,$1":"$2 }}' ../2_boundaries_with_DI/combined_boundary.uniq.txt ) \
      -b $file \
  -wo > ${sample}.DI.overlap.txt
  done

Rscript calc_DI_delta.r

Rscript define_boundary_by_DI_delta.r
Rscript find_dynamic_boundary_DI_delta.r

#Dynamic txt files to bed with 0, 25kb, 50kb, 100kb extension on both sides of boundary
cd 4_dynamic_boundary/
for file in `ls *txt` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' $file > ${xbase/.txt/.bed} ; done
for file in `ls *txt` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '{print $1,$2-25000,$2+25000,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' $file > ${xbase/.txt/_25kb.bed} ; done
for file in `ls *txt` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '{print $1,$2-50000,$2+50000,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' $file > ${xbase/.txt/_50kb.bed} ; done
for file in `ls *txt` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '{print $1,$2-100000,$2+100000,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' $file > ${xbase/.txt/_100kb.bed} ; done

#TAD size/number and comparison with domain caller
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000}' /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/2_TAD/Bing/diff_bound/2_boundaries_with_DI/boundary.all.uniq.txt > domain_caller.boundary.all.uniq.bed
awk -v FS="\t" -v OFS="\t" '{print $1,$2-50000,$2+50000}' ../2_boundaries_with_DI/boundary.all.uniq.txt > 3DNetMod.boundary.all.uniq.bed
tail -n +2 ../3_non_redundant_boundary/GSC_combined_boundary.DI_cutoff.100.all.tsv | sed 's/"//g' | awk -v FS="\t" -v OFS="\t" '{print $1,$2-50000,$2+50000}' > GSC_combined_boundary.DI_cutoff.100.all.bed
tail -n +2 ../3_non_redundant_boundary/GSC_combined_boundary.DI_cutoff.25.all.tsv | sed 's/"//g' | awk -v FS="\t" -v OFS="\t" '{print $1,$2-50000,$2+50000}' > GSC_combined_boundary.DI_cutoff.25.all.bed
tail -n +2 ../3_non_redundant_boundary/GSC_combined_boundary.DI_cutoff.0.all.tsv | sed 's/"//g' | awk -v FS="\t" -v OFS="\t" '{print $1,$2-50000,$2+50000}' > GSC_combined_boundary.DI_cutoff.0.all.bed

wc -l domain_caller.boundary.all.uniq.bed
148975 domain_caller.boundary.all.uniq.bed
bedtools intersect -a domain_caller.boundary.all.uniq.bed -b GSC_combined_boundary.DI_cutoff.100.all.bed -u | wc -l
102821
bedtools intersect -a domain_caller.boundary.all.uniq.bed -b GSC_combined_boundary.DI_cutoff.25.all.bed -u | wc -l
121262
bedtools intersect -a domain_caller.boundary.all.uniq.bed -b GSC_combined_boundary.DI_cutoff.0.all.bed -u | wc -l

wc -l GSC_combined_boundary.DI_cutoff.25.all.bed
11082 GSC_combined_boundary.DI_cutoff.25.all.bed
wc -l GSC_combined_boundary.DI_cutoff.100.all.bed
6222 GSC_combined_boundary.DI_cutoff.100.all.bed
bedtools intersect -b domain_caller.boundary.all.uniq.bed -a GSC_combined_boundary.DI_cutoff.100.all.bed -u | wc -l
4613
bedtools intersect -b domain_caller.boundary.all.uniq.bed -a GSC_combined_boundary.DI_cutoff.25.all.bed -u | wc -l
6832

awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000,$1,$3-10000,$3,$4}' ../2_boundaries_with_DI/combined_tads.raw.sorted.txt | bedtools pairtobed -a stdin -b GSC_combined_boundary.DI_cutoff.100.all.bed -type both | cut -f1,2,6,7 | sort -k1,1 -k2,2n -k3,3n | uniq > combined_tads.DI_cutoff.100.sorted.txt
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000,$1,$3-10000,$3,$4}' ../2_boundaries_with_DI/combined_tads.raw.sorted.txt | bedtools pairtobed -a stdin -b GSC_combined_boundary.DI_cutoff.25.all.bed -type both | cut -f1,2,6,7 | sort -k1,1 -k2,2n -k3,3n | uniq > combined_tads.DI_cutoff.25.sorted.txt


#Get TADs with two 'valid' boundaries
cd ..
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000,$1,$3-10000,$3}' 2_boundaries_with_DI/combined_tads.raw.sorted.txt | bedtools pairtobed -a stdin -b 4_dynamic_boundary/all_boundary_withDI_new.bed -type both | cut -f1,2,6 | sort -k1,1 -k2,2n -k3,3n  > 4_dynamic_boundary/all_TADs.bed
#R notebook for all valid TADs
cd 4_dynamic_boundary/
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000,$1,$3-10000,$3}' all_TADs_merged.bed | bedtools pairtobed -a stdin -b dynamic_sig_withDI_new_25kb.bed -type either | cut -f1,2,6 | sort -k1,1 -k2,2n -k3,3n > dynamic_TADs.bed
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000,$1,$3-10000,$3}' all_TADs_merged.bed | bedtools pairtobed -a stdin -b dynamic_strengthened_sig_withDI_new_25kb.bed -type either | cut -f1,2,6 | sort -k1,1 -k2,2n -k3,3n > strengthened_TADs.bed
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000,$1,$3-10000,$3}' all_TADs_merged.bed | bedtools pairtobed -a stdin -b dynamic_weakened_sig_withDI_new_25kb.bed -type either | cut -f1,2,6 | sort -k1,1 -k2,2n -k3,3n > weakened_TADs.bed
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$2+10000,$1,$3-10000,$3}' all_TADs_merged.bed | bedtools pairtobed -a stdin -b dynamic_sig_withDI_new_25kb.bed -type neither | cut -f1,2,6 | sort -k1,1 -k2,2n -k3,3n > not_dynamic_TADs.bed
#Overlap with RNA
for file in `ls *DEGs.bed` ; do xbase=$(basename $file) ; tail -n +2 $file | bedtools intersect -a ../4_dynamic_boundary/strengthened_TADs.bed -b stdin -wo > ${xbase/.bed/_stronger_TADs.bed} ; done
for file in `ls *DEGs.bed` ; do xbase=$(basename $file) ; tail -n +2 $file | bedtools intersect -a ../4_dynamic_boundary/not_dynamic_TADs.bed -b stdin -wo > ${xbase/.bed/_notdynamic_TADs.bed} ; done
for file in `ls *DEGs.bed` ; do xbase=$(basename $file) ; tail -n +2 $file | bedtools intersect -a ../4_dynamic_boundary/weakened_TADs.bed -b stdin -wo > ${xbase/.bed/_weaker_TADs.bed} ; done

#5_aggregate_boundaries
for bed in `ls *bed` ; do for file in `ls /scratch/mmoore/Epitherapy3D/data/HiC/condition_inter30_hic/*hic` ; do xbase=$(basename $file) ; xbase1=$(basename $bed) ; xbase2=${xbase1/.bed/}; echo fanc aggregate ${file}@10kb@VC_SQRT $bed ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg} -p ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg.png} -m ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg.txt} -a 200000 -e -l --vmin=-0.5 --vmax=0.5 ; done ; done > fanc_aggregate_commands.txt
parallel_GNU -j 5 < fanc_aggregate_commands.txt

###WGBS
for file in `ls ../6_track_mergedCG/*methylC.gz` ; do xbase=$(basename $file) ; zcat $file | awk -v FS="\t" -v OFS="\t" '$7 > 4 {print $1,$2,$3,$5}' > ${xbase/.methylC.gz/.bedgraph} ; done
for file in `ls *bedgraph` ; do xbase=$(basename $file) ; echo bedGraphToBigWig $file /scratch/devtools/mmoore/genomes/private/hg38_all_chromosomes.size ${xbase/.bedgraph/.bw} ; done > bg_to_bw.txt
parallel_GNU -j 10 < bg_to_bw.txt

#computeMatrix scale-regions -S ./*bw -R ../4_dynamic_boundary/not_dynamic_sig_withDI_new_100kb.bed ../4_dynamic_boundary/dynamic_weakened_sig_withDI_new_100kb.bed  ../4_dynamic_boundary/dynamic_strengthened_sig_withDI_new_100kb.bed -o mCG_TAD.matrix.gz --outFileNameMatrix mCG_TAD.matrix.tab -p 10
computeMatrix reference-point -S ./*.bw -R ../4_dynamic_boundary/not_dynamic_sig_withDI_new_100kb.bed ../4_dynamic_boundary/dynamic_weakened_sig_withDI_new_100kb.bed  ../4_dynamic_boundary/dynamic_strengthened_sig_withDI_new_100kb.bed -o mCG_TAD_new.matrix.gz --outFileNameMatrix mCG_TAD_new.matrix.tab -p 20 --referencePoint "center" -a 100000 -b 100000 -bs 5000 --smartLabels
#plotHeatmap -m mCG_TAD.matrix.gz -o mCG_TAD.png 
plotHeatmap -m mCG_TAD_new.matrix.gz -o mCG_TAD_new.png --xAxisLabel "" --refPointLabel "Boundary" --regionsLabel "Stable Boundaries" "DMSO Boundaries" "DP Boundaries" --colorList 'blue, white, red' --yMin 0 --yMax 100

###CTCF
for file in `ls ../1_diffpeaks/*bed` ; do xbase=$(basename $file) ; bedtools intersect -a peaks.summits -b $file -u > ${xbase/.bed/_summits.bed} ; done
for file in `ls ../1_diffpeaks/*bed` ; do xbase=$(basename $file) ; bedtools intersect -a peaks.summits -b $file -u | bedtools intersect -a stdin -b /scratch/devtools/mmoore/genomes/repeatmasker/UCSC_rmsk.bed -u > ${xbase/.bed/_summits_TE.bed} ; done

echo -e "peaks\tTE\tall" > GSC_summits_TE_counts.txt
for file in `ls *summits.bed` ; do xbase=$(basename $file) ; name=${xbase/_peaks_summits.bed/} ; a1=$(cat $file | wc -l) a2=$(cat ${file/.bed/_TE.bed} | wc -l) ; echo -e "${name}\t${a2}\t${a1}" >> GSC_summits_TE_counts.txt ; done

echo -e "peaks\tTE\tall" > GSC_CTCF_TE_counts.txt
for file in `ls ../1_diffpeaks/*peaks.bed` ; do xbase=$(basename $file) ; name=${xbase/.bed/_CTCF_TE} ; a1=$(cat $file | wc -l) a2=$(intersectBed -a /scratch/devtools/mmoore/genomes/repeatmasker/UCSC_rmsk.bed -b ../1_diffpeaks/CTCF_hg38_all.bed -u | intersectBed -a $file -b stdin -u | wc -l) ; echo -e "${name}\t${a2}\t${a1}" >> GSC_CTCF_TE_counts.txt ; done

for file in `ls ../6_TAD_ctcf/*25kb.bed` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3}' $file | bedtools coverage -a stdin -b ../1_diffpeaks/TE_CTCF_hg38.bed > ${xbase/.bed/_CTCF_TE_coverage.bed} ; done

for file in `ls ../1_diffpeaks/*peaks.bed` ; do xbase=$(basename $file) ; intersectBed -a /scratch/devtools/mmoore/genomes/repeatmasker/UCSC_rmsk.bed -b ../1_diffpeaks/CTCF_hg38_all.bed -u | intersectBed -b $file -a stdin -u > ${xbase/.bed/_CTCF_TEs.bed} ; done

#for file in `ls ../1_diffpeaks/*bed` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3}' $file | bedtools coverage -a stdin -b /scratch/devtools/mmoore/genomes/repeatmasker/UCSC_rmsk.bed > ${xbase/.bed/_coverage.bed} ; done

for file in `ls ../1_diffpeaks/*peaks.bed` ; do for file1 in /scratch/mmoore/Epitherapy3D/data/wgbs/methylcall/6_track_mergedCG/*CG.methylC.gz ; do xbase=$(basename $file) ; xbase1=$(basename $file1) ; name=${xbase1/Trimmed_WGBS_/} ; awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3}' $file | bedtools intersect -a stdin -b $file1 -wo > ${xbase/.bed/}_${name/.gz/} ; done ; done
for file in `ls ../1_diffpeaks/*peaks.bed` ; do for file1 in /scratch/mmoore/Epitherapy3D/data/wgbs/methylcall/6_track_mergedCG/*CG.methylC.gz ; do xbase=$(basename $file) ; xbase1=$(basename $file1) ; name=${xbase1/Trimmed_WGBS_/} ; awk -v FS="\t" -v OFS="\t" '{print $1,$2-1000,$3+1000}' $file | bedtools intersect -a stdin -b $file1 -wo > ${xbase/.bed/}_${name/.CG.methylC.gz/_1kb.CG.methylC} ; done ; done
for file in `ls ../1_diffpeaks/*peaks.bed` ; do for file1 in /scratch/mmoore/Epitherapy3D/data/wgbs/methylcall/6_track_mergedCG/*CG.methylC.gz ; do xbase=$(basename $file) ; xbase1=$(basename $file1) ; name=${xbase1/Trimmed_WGBS_/} ; awk -v FS="\t" -v OFS="\t" '{print $1,$2-2000,$3+2000}' $file | bedtools intersect -a stdin -b $file1 -wo > ${xbase/.bed/}_${name/.CG.methylC.gz/_2kb.CG.methylC} ; done ; done


computeMatrix reference-point -S ./*bw -R ../1_diffpeaks/GSC_CTCF_nonSig_peaks.bed ../1_diffpeaks/GSC_CTCF_downSig_peaks.bed ../1_diffpeaks/GSC_CTCF_upSig_peaks.bed -o mCG_CTCF.matrix.gz --outFileNameMatrix mCG_CTCF.matrix.tab -p 10 -bs 1500 --referencePoint "center" -a 6000 -b 6000
plotHeatmap -m mCG_CTCF.matrix.gz -o mCG_CTCF.pdf --regionsLabel "Stable Peaks" "Lost Peaks" "Gained Peaks" --xAxisLabel "" --refPointLabel "center" --colorMap 'bwr' --yMin 0 --yMax 100

computeMatrix reference-point -S ./*bw -R ../1_diffpeaks/GSC_CTCF_nonSig_peaks.bed ../1_diffpeaks/GSC_CTCF_downSig_peaks.bed ../1_diffpeaks/GSC_CTCF_upSig_peaks.bed -o mCG_CTCF_new.matrix.gz --outFileNameMatrix mCG_CTCF_new.matrix.tab -bs 1500 --referencePoint "center" -a 6000 -b 6000 -p 20 --minThreshold 0
plotHeatmap -m mCG_CTCF_new.matrix.gz -o mCG_CTCF_new.pdf --regionsLabel "Stable Peaks" "Lost Peaks" "Gained Peaks" --xAxisLabel "" --refPointLabel "center" --colorMap 'bwr' --yMin 0 --yMax 100

echo -e "bound_type\tCTCFdown\tCTCFdownONLY\tCTCFup\tCTCFupONLY\tCTCFboth\tCTCFstable\tnoCTCF\tTotal" > ctcf_TAD_change_overlap.txt
for file in $(ls ./*25kb.bed); do
  name=$(basename $file)
  bound=${name%.bed}
  a1=$(intersectBed -a $file -b ../2_allpeaks_anno/GSC_CTCF_downSig_peaks.bed -u |wc -l)
  a2=$(intersectBed -a $file -b ../2_allpeaks_anno/GSC_CTCF_downSig_peaks.bed -u | bedtools intersect -a stdin -b ../2_allpeaks_anno/GSC_CTCF_upSig_peaks.bed -v | wc -l)
  a3=$(intersectBed -a $file -b ../2_allpeaks_anno/GSC_CTCF_upSig_peaks.bed -u |wc -l)
  a4=$(intersectBed -a $file -b ../2_allpeaks_anno/GSC_CTCF_upSig_peaks.bed -u | bedtools intersect -a stdin -b ../2_allpeaks_anno/GSC_CTCF_downSig_peaks.bed -v | wc -l)
  a5=$(intersectBed -a $file -b ../2_allpeaks_anno/GSC_CTCF_upSig_peaks.bed -u | bedtools intersect -a stdin -b ../2_allpeaks_anno/GSC_CTCF_downSig_peaks.bed -u | wc -l)
  a6=$(intersectBed -a $file -b ../2_allpeaks_anno/GSC_CTCF_nonSig_peaks.bed -u | bedtools intersect -a stdin -b ../2_allpeaks_anno/GSC_CTCF_downSig_peaks.bed -v | bedtools intersect -a stdin -b ../2_allpeaks_anno/GSC_CTCF_upSig_peaks.bed -v | wc -l)
  a7=$(intersectBed -a $file -b ../2_allpeaks_anno/GSC_CTCF_nonSig_peaks.bed -v |  bedtools intersect -a stdin -b ../2_allpeaks_anno/GSC_CTCF_downSig_peaks.bed -v | bedtools intersect -a stdin -b ../2_allpeaks_anno/GSC_CTCF_upSig_peaks.bed -v | wc -l)
  a8=$(cat $file | wc -l)
  echo -e "$bound\t${a1}\t${a2}\t${a3}\t${a4}\t${a5}\t${a6}\t${a7}\t${a8}" |tee -a ctcf_TAD_change_overlap.txt
  done

for file in $(ls ./*25kb.bed); do
  name=$(basename $file)
  bound=${name%_with*}
  intersectBed -a $file -b ../1_diffpeaks/GSC_CTCF_nonSig_peaks.bed -wo | bedtools intersect -a stdin -b ../1_diffpeaks/GSC_CTCF_downSig_peaks.bed -v | 
  bedtools intersect -a stdin -b ../1_diffpeaks/GSC_CTCF_upSig_peaks.bed -v > ${bound}_CTCF_stable.bed
  done

###Loops
cat ../../3_CTCF/1_diffpeaks/GSC_CTCF_downSig_peaks.bed ../../3_CTCF/1_diffpeaks/GSC_CTCF_upSig_peaks.bed > all_diffpeaks.bed
for file in *diffloops.bedpe ; do xbase=$(basename $file) ; bedtools pairtobed -a $file -b all_diffpeaks.bed > ${xbase/.bedpe/_ctcf.bedpe} ; done
for file in *diffloops.bedpe ; do xbase=$(basename $file) ; bedtools pairtobed -a $file -b all_diffpeaks.bed | cut -f1-6 | sort -k1,1 -k2,2n | uniq > ${xbase/.bedpe/_ctcf_uniq.bedpe} ; done
for file in *diffloops.bedpe ; do xbase=$(basename $file) ; cell=${xbase/_diffloops.bedpe/} ; tail -n +2 ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | bedtools pairtobed -a $file -b stdin > ${xbase/.bedpe/_genes.bedpe} ; done
for file in *diffloops.bedpe ; do xbase=$(basename $file) ; cell=${xbase/_diffloops.bedpe/} ; tail -n +2 ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | bedtools pairtobed -a $file -b stdin | cut -f1-6 | sort -k1,1 -k2,2n | uniq > ${xbase/.bedpe/_genes_uniq.bedpe} ; done
for file in *diffloops.bedpe ; do xbase=$(basename $file) ; cell=${xbase/_diffloops.bedpe/} ; awk '$7 > 1 || $7 < -1 && $10 < 0.05 {print $0}' ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | tail -n +2 |bedtools pairtobed -a $file -b stdin > ${xbase/.bedpe/_deg.bedpe} ; done
for file in *diffloops.bedpe ; do xbase=$(basename $file) ; cell=${xbase/_diffloops.bedpe/} ; awk '$7 > 1 || $7 < -1 && $10 < 0.05 {print $0}' ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | tail -n +2 |bedtools pairtobed -a $file -b stdin | cut -f1-6 | sort -k1,1 -k2,2n | uniq > ${xbase/.bedpe/_deg_uniq.bedpe} ; done
for file in *diffloops_ctcf_uniq.bedpe ; do xbase=$(basename $file) ; cell=${xbase/_diffloops_ctcf_uniq.bedpe/} ; tail -n +2 ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | bedtools pairtobed -a $file -b stdin > ${xbase/.bedpe/_genes.bedpe} ; done
for file in *diffloops_ctcf_uniq.bedpe ; do xbase=$(basename $file) ; cell=${xbase/_diffloops_ctcf_uniq.bedpe/} ; tail -n +2 ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | bedtools pairtobed -a $file -b stdin | cut -f1-6 | sort -k1,1 -k2,2n | uniq > ${xbase/.bedpe/_genes_uniq.bedpe} ; done
for file in *diffloops_ctcf_uniq.bedpe ; do xbase=$(basename $file) ; cell=${xbase/_diffloops_ctcf_uniq.bedpe/} ; awk '$7 > 1 || $7 < -1 && $10 < 0.05 {print $0}' ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | tail -n +2 |bedtools pairtobed -a $file -b stdin > ${xbase/.bedpe/_deg.bedpe} ; done
for file in *diffloops_ctcf_uniq.bedpe ; do xbase=$(basename $file) ; cell=${xbase/_diffloops_ctcf_uniq.bedpe/} ; awk '$7 > 1 || $7 < -1 && $10 < 0.05 {print $0}' ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | tail -n +2 |bedtools pairtobed -a $file -b stdin | cut -f1-6 | sort -k1,1 -k2,2n | uniq > ${xbase/.bedpe/_deg_uniq.bedpe} ; done

echo -e "cell\tCTCF\tCTCFONLY\tDEG\tDEGONLY\tCTCF_DEG\tneither\tTotal" > ctcf_loop_change_overlap.txt
for file in $(ls *diffloops.bedpe); do
  xbase=$(basename $file)
  cell=${xbase/_diffloops.bedpe/}
  a1=$(cat ${cell}_diffloops_ctcf_uniq.bedpe |wc -l)
  a2=$(awk '$7 > 1 || $7 < -1 && $10 < 0.05 {print $0}' ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | tail -n +2 | bedtools pairtobed -a ${cell}_diffloops_ctcf_uniq.bedpe -b stdin -type neither | cut -f1-6 | sort -k1,1 -k2,2n | uniq | wc -l)
  a3=$(cat ${cell}_diffloops_deg_uniq.bedpe |wc -l)
  a4=$(bedtools pairtobed -a ${cell}_diffloops_deg_uniq.bedpe -b all_diffpeaks.bed -type neither | cut -f1-6 | sort -k1,1 -k2,2n | uniq | wc -l)
  a5=$(cat ${cell}_diffloops_ctcf_uniq_deg_uniq.bedpe | wc -l)
  a6=$(awk '$7 > 1 || $7 < -1 && $10 < 0.05 {print $0}' ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | tail -n +2 | bedtools pairtobed -a $file -b stdin -type neither | cut -f1-6 | sort -k1,1 -k2,2n | uniq | bedtools pairtobed -a stdin -b all_diffpeaks.bed -type neither | cut -f1-6 | sort -k1,1 -k2,2n | uniq | wc -l)
  a7=$(cat $file | wc -l)
  echo -e "$cell\t${a1}\t${a2}\t${a3}\t${a4}\t${a5}\t${a6}\t${a7}" |tee -a ctcf_loop_change_overlap.txt
  done

for file in `ls ../diffloops_old/new_unsampled/*dloops.bedpe` ; do xbase=$(basename $file) ; bedtools pairtobed -a $file -b ../../2_TAD/3DNetMod/diff_bound/4_dynamic_boundary/dynamic_sig_withDI_new.bed -type ispan | cut -f1-6 | sort -k1,1 -k2,2n | uniq > ${xbase/.bedpe/_span_diffbound.bedpe} ; done

echo -e "cell\tDMSO_cross\tDP_cross\tAll_cross\tDMSO_Total\tDP_Total\tAll_Total" > bound_loop_change_overlap.txt
for file in $(ls ../diffloops_old/new_unsampled/*DMSOvsDP_dloops.bedpe); do
  xbase=$(basename $file)
  cell=${xbase/_DMSOvsDP_dloops.bedpe/}
  a1=$(cat ${cell}_DMSOvsDP_${cell}_DMSO_specific_dloops_span_diffbound.bedpe |wc -l)
  a2=$(cat ${cell}_DMSOvsDP_${cell}_DP_specific_dloops_span_diffbound.bedpe |wc -l)
  a3=$(cat ${cell}_DMSOvsDP_dloops_span_diffbound.bedpe |wc -l)
  a4=$(cat ../diffloops_old/new_unsampled/${cell}_DMSOvsDP_${cell}_DMSO_specific_dloops.bedpe |wc -l)
  a5=$(cat ../diffloops_old/new_unsampled/${cell}_DMSOvsDP_${cell}_DP_specific_dloops.bedpe |wc -l)
  a6=$(cat $file | wc -l)
  echo -e "$cell\t${a1}\t${a2}\t${a3}\t${a4}\t${a5}\t${a6}" |tee -a bound_loop_change_overlap.txt
  done

for file in $(ls ../diffloops_old/new_unsampled/*_specific_dloops.bedpe); do
  xbase=$(basename $file)
  cell=${xbase/_DMSOvsDP_*/}
  tail -n +2 ../../rnaseq/1_DESeq2/${cell}_DPvsDMSO_shrinkage_DEGs.bed | bedtools pairtobed -a $file -b stdin > ${xbase/.bedpe/_genes.bedpe}
done

for file in `ls ../4_loop_sep_deg/*` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '$13 > 1 || $13 < -1 && $16<0.05 {print $1,$2,$3,$4,$5,$6}' $file | sort -k1,1 -k2,2n | uniq | wc -l ; done
wc -l ../diffloops_old/new_unsampled/*specific*bedpe

for file in `ls ../diffloops_old/new_unsampled/*_specific_dloops.bedpe` ; do xbase=$(basename $file) ; cat $file | sed 's/chr//g' > ${xbase/.bedpe/_nochr.bedpe} ; done
for bedpe in `ls *_nochr.bedpe` ; do for file in `ls /scratch/mmoore/Epitherapy3D/data/HiC/condition_inter30_hic/*hic` ; do xbase=$(basename $file) ; xbase1=$(basename $bedpe) ; xbase2=${xbase1/_dloops.bedpe/} ; echo fanc aggregate ${file}@10kb@VC_SQRT $bedpe ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg} -p ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg.png} -m ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg.txt} --loops --vmin=-1 --vmax=1 ; done ; done > fanc_aggregate_commands.txt
parallel_GNU -j 3 < fanc_aggregate_commands.txt

cd ../6_shared_diff/
for bedpe in `ls *bedpe` ; do for file in `ls /scratch/mmoore/Epitherapy3D/data/HiC/condition_inter30_hic/*hic` ; do xbase=$(basename $file) ; xbase1=$(basename $bedpe) ; xbase2=${xbase1/.bedpe/} ; echo fanc aggregate ${file}@10kb@VC_SQRT $bedpe ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg} -p ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg.png} -m ${xbase2}_${xbase/.hic/_10kb_VC_SQRT.agg.txt} --loops --vmin=-2 --vmax=2 ; done ; done > fanc_aggregate_commands.txt
parallel_GNU -j 3 < fanc_aggregate_commands.txt

cd ../7_merged_loop/
echo -e "loopId\tchrA\tstartA\tendA\tchrB\tstartB\tendB\tdistance(bp)" > GSC_shared_differential_loops.txt
awk -v FS="\t" -v OFS="\t" '{print "loop_chr"$1"-chr"$4"-"NR,"chr"$1,$2,$3,"chr"$4,$5,$6,(($5+$6)/2-($2+$3)/2)}' GSC_shared_differential_loops.bedpe | perl -lane '$F[7]=sprintf("%.15f",$F[7]) if $F[7]=~/e/i; $F[7]=~s/\.?0+$//; print join("\t",@F)' >> GSC_shared_differential_loops.txt

cLoops2 quant -d ../../../4_loop/cloops/B36_DMSO -loops GSC_shared_differential_loops.txt -o B36_DMSO_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/B36_DP -loops GSC_shared_differential_loops.txt -o B36_DP_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/B49_DMSO -loops GSC_shared_differential_loops.txt -o B49_DMSO_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/B49_DP -loops GSC_shared_differential_loops.txt -o B49_DP_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/B66_DMSO -loops GSC_shared_differential_loops.txt -o B66_DMSO_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/B66_DP -loops GSC_shared_differential_loops.txt -o B66_DP_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/NHA_DMSO -loops GSC_shared_differential_loops.txt -o NHA_DMSO_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/NHA_DP -loops GSC_shared_differential_loops.txt -o NHA_DP_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/qNHA_DMSO -loops GSC_shared_differential_loops.txt -o qNHA_DMSO_diffloops_quant -p 25
cLoops2 quant -d ../../../4_loop/cloops/qNHA_DP -loops GSC_shared_differential_loops.txt -o qNHA_DP_diffloops_quant -p 25

for file in `ls *loops.txt` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '{print $2,$3,$4,$5,$6,$7,$1,$16,$21}' $file > ${xbase/txt/bedpe} ; done

#TE enrichment script

tail -n +2 GSC_shared_differential_loops.bedpe | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3"\n"$4,$5,$6}' > GSC_shared_differential_anchors.bed

zcat /scratch/yliang/HNSCC/analysis/2_CPTAC_TSTEA_candidates/3_test3_real_TCGA_data/teprof3_v3.1.8_gencodeV26_curated/cell_line_GBM/8_tumor_specificity_summary/teprof3_output_quantification.TE.cell_line_GBM.filtered.tsv.gz | awk -v FS="\t" -v OFS="\t" '$2 >=1 && $9>=1 && $9>=2*$10 {print $0}' > teprof3_output_quantification.TE.cell_line_GBM.filtered.tsv
for i in "B36" "B66" "B49" ; do
  bedtools pairtobed -a ../diffloops_old/new_unsampled/${i}_DMSOvsDP_${i}_DMSO_specific_dloops.bedpe -b teprof3_tumor_specific_DP_specific_TSTET.bed | grep $i >> diff_loop_indiv_teprof3_tumor_specific_DP_specific_TSTET.txt
  bedtools pairtobed -a ../diffloops_old/new_unsampled/${i}_DMSOvsDP_${i}_DP_specific_dloops.bedpe -b teprof3_tumor_specific_DP_specific_TSTET.bed | grep $i >> diff_loop_indiv_teprof3_tumor_specific_DP_specific_TSTET.txt
done

echo -e "cell\tDMSO_loop_TET\tDP_loop_TET\tTET" > loop_tet_overlap_indiv.txt
for i in "B36" "B66" "B49" ; do
  a1=$(bedtools pairtobed -a ../diffloops_old/new_unsampled/${i}_DMSOvsDP_${i}_DMSO_specific_dloops.bedpe -b teprof3_tumor_specific_DP_specific_TSTET.bed | grep $i | cut -f10 | sort | uniq | wc -l)
  a2=$(bedtools pairtobed -a ../diffloops_old/new_unsampled/${i}_DMSOvsDP_${i}_DP_specific_dloops.bedpe -b teprof3_tumor_specific_DP_specific_TSTET.bed | grep $i | cut -f10 | sort | uniq | wc -l)
  a3=$(cat teprof3_tumor_specific_DP_specific_TSTET.bed | grep $i | cut -f4 | sort | uniq | wc -l)
  echo -e "$i\t${a1}\t${a2}\t${a3}" |tee -a loop_tet_overlap_indiv.txt
  done

#browser
cat ../../../cutrun/1_diffbind/GSC_CTCF_downSig_peaks.bed ../../../cutrun//1_diffbind/GSC_CTCF_upSig_peaks.bed | sort -k1,1 -k2,2n > GSC_CTCF_diffSig_peaks.bed
bigwigAverage -b /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B36-DMSO1-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B36-DMSO2-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw -o B36_DMSO_CTCF.bw -p 15
bigwigAverage -b /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B36-DP1-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B36-DP2-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw -o B36_DP_CTCF.bw -p 15
bigwigAverage -b /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B49-DMSO1-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B49-DMSO2-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw -o B49_DMSO_CTCF.bw -p 15
bigwigAverage -b /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B49-DP1-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B49-DP2-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw -o B49_DP_CTCF.bw -p 15
bigwigAverage -b /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B66-DMSO1-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B66-DMSO2-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw -o B66_DMSO_CTCF.bw -p 15
bigwigAverage -b /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B66-DP1-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw /scratch/mmoore/Epitherapy3D/data/cutrun_all_samples/bigwig_normalized/RLE/Trimmed_WangT_B66-DP2-CTCF-cutrun_R1.fastq_nochrM_nodup_qffilter.bw -o B66_DP_CTCF.bw -p 15

computeMatrix reference-point -S ../0_norm_bw/*bw -R ../1_diffbind/GSC_CTCF_downSig_peaks.bed ../1_diffbind/GSC_CTCF_upSig_peaks.bed -a 3000 -b 3000 --referencePoint center -p 15 -o GSC_CTCF_peaks.matrix.gz --missingDataAsZero
plotHeatmap -m GSC_CTCF_peaks.matrix.gz -out Fig_3A.pdf --colorMap Greens Greens Greens Greens Greens Greens --whatToShow 'heatmap and colorbar' --regionsLabel "Lost CTCF Sites (n=3893)" "Gained CTCF Sites (n=4785)" --samplesLabel "B36 DMSO" "B36 DP" "B49 DMSO" "B49 DP" "B66 DMSO" "B66 DP" --refPointLabel "center" --xAxisLabel "" --heatmapHeight 10 --heatmapWidth 6

computeMatrix reference-point -S ../0_norm_bw/*bw -R ../1_diffbind/GSC_CTCF_nonSig_peaks.bed ../1_diffbind/GSC_CTCF_downSig_peaks.bed ../1_diffbind/GSC_CTCF_upSig_peaks.bed -a 3000 -b 3000 --referencePoint center -p 15 -o GSC_CTCF_peaks_with_nonsig.matrix.gz --missingDataAsZero
plotProfile -m GSC_CTCF_peaks_with_nonsig.matrix.gz -out Fig_3B.pdf --regionsLabel "Unchanged CTCF Sites (n=43,287)" "Lost CTCF Sites (n=3893)" "Gained CTCF Sites (n=4785)" --samplesLabel "B36 DMSO" "B36 DP" "B49 DMSO" "B49 DP" "B66 DMSO" "B66 DP" --refPointLabel "center" --perGroup --colors "#33A02C" "#B2DF8A" "#FF7F00" "#FDBF6F" "#1F78B4" "#A6CEE3" --numPlotsPerRow 2 

for file in `ls ../1_diffbind/*peaks.bed` ; do xbase=$(basename $file) ; awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"_"$3,".","+"}' $file > ${xbase/.bed/_homer.bed} ; done



###snm3c-seq
computeMatrix reference-point -S B36*bw -R /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/2_TAD/3DNetMod/diff_bound/4_dynamic_boundary/dynamic_weakened_sig_withDI_new.bed /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/2_TAD/3DNetMod/diff_bound/4_dynamic_boundary/dynamic_strengthened_sig_withDI_new.bed -a 200000 -b 200000 --referencePoint center -p 20 -o B36_dynamic_split_matrix_new.gz --missingDataAsZero --skipZeros
plotHeatmap -m B36_dynamic_split_matrix_new.gz -out B36_dynamic_split_matrix_new.pdf --regionsLabel "DMSO Boundaries" "DP Boundaries" --samplesLabel "B36_0" "B36_1" "B36_2" "B36_3" "B36_4" "B36_5" "B36_6" --refPointLabel "Boundary" --whatToShow 'heatmap and colorbar' --colorMap 'bwr' --xAxisLabel ""
plotProfile -m B36_dynamic_split_matrix_new.gz -out B36_dynamic_split_matrix_new_profile.pdf --regionsLabel "DMSO Boundaries" "DP Boundaries" --samplesLabel "B36_0" "B36_1" "B36_2" "B36_3" "B36_4" "B36_5" "B36_6" --refPointLabel "Boundary" --outFileNameData B36_dynamic_split_matrix_new_profile.tab
computeMatrix reference-point -S B49*bw -R /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/2_TAD/3DNetMod/diff_bound/4_dynamic_boundary/dynamic_weakened_sig_withDI_new.bed /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/2_TAD/3DNetMod/diff_bound/4_dynamic_boundary/dynamic_strengthened_sig_withDI_new.bed -a 200000 -b 200000 --referencePoint center -p 20 -o B49_dynamic_split_matrix_new.gz --missingDataAsZero --skipZeros
plotHeatmap -m B49_dynamic_split_matrix_new.gz -out B49_dynamic_split_matrix_new.pdf --regionsLabel "DMSO Boundaries" "DP Boundaries" --samplesLabel "B49_0" "B49_1" "B49_2" "B49_3" "B49_4" "B49_5" "B49_6" --refPointLabel "Boundary" --whatToShow 'heatmap and colorbar' --colorMap 'bwr' --xAxisLabel ""
plotProfile -m B49_dynamic_split_matrix_new.gz -out B49_dynamic_split_matrix_new_profile.pdf --regionsLabel "DMSO Boundaries" "DP Boundaries" --samplesLabel "B49_0" "B49_1" "B49_2" "B49_3" "B49_4" "B49_5" "B49_6" --refPointLabel "Boundary" --outFileNameData B49_dynamic_split_matrix_new_profile.tab
computeMatrix reference-point -S B66*bw -R /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/2_TAD/3DNetMod/diff_bound/4_dynamic_boundary/dynamic_weakened_sig_withDI_new.bed /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/2_TAD/3DNetMod/diff_bound/4_dynamic_boundary/dynamic_strengthened_sig_withDI_new.bed -a 200000 -b 200000 --referencePoint center -p 20 -o B66_dynamic_split_matrix_new.gz --missingDataAsZero --skipZeros
plotHeatmap -m B66_dynamic_split_matrix_new.gz -out B66_dynamic_split_matrix_new.pdf --regionsLabel "DMSO Boundaries" "DP Boundaries" --samplesLabel "B66_0" "B66_1" "B66_2" "B66_3" "B66_4" "B66_5" "B66_6" --refPointLabel "Boundary" --whatToShow 'heatmap and colorbar' --colorMap 'bwr' --xAxisLabel ""
plotProfile -m B66_dynamic_split_matrix_new.gz -out B66_dynamic_split_matrix_new_profile.pdf --regionsLabel "DMSO Boundaries" "DP Boundaries" --samplesLabel "B66_0" "B66_1" "B66_2" "B66_3" "B66_4" "B66_5" "B66_6" --refPointLabel "Boundary" --outFileNameData B66_dynamic_split_matrix_new_profile.tab


for file in `ls ../3_CTCF/1_diffpeaks/*peaks.bed` ; do xbase=$(basename $file) ; xbase1=${xbase/Sig_peaks.bed/} ; xbase2=${xbase1/GSC_CTCF_/} ; awk -v FS="\t" -v OFS="\t" -v type=$xbase2 '{print $1,$2,$3,type"_"NR}' $file >> all_CTCF_anno.bed ; done
cat all_CTCF_anno.bed | sort -k1,1 -k2,2n > all_CTCF_anno_sorted.bed
awk -v FS="\t" -v OFS="\t" '{print $1,$2-1000,$3+1000,$4}' all_CTCF_anno_sorted.bed > all_CTCF_1k_anno_sorted.bed
awk -v FS="\t" -v OFS="\t" '{print $1,$2-500,$3+500,$4}' all_CTCF_anno_sorted.bed > all_CTCF_500bp_anno_sorted.bed
awk -v FS="\t" -v OFS="\t" '{print $1, ($3-$2<1000 ? $2-500 : $2), ($3-$2<1000 ? $3+500 : $3), $4}' all_CTCF_anno_sorted.bed > all_CTCF_scaled_anno_sorted.bed

allcools generate-dataset  \
--allc_table AllcPaths.tsv \
--output_path GSC_epitherapy_2kb_new.mcds \
--chrom_size_path /scratch/devtools/mmoore/genomes/snm3c/hg38/chrom_sizes.txt \
--obs_dim cell  \
--cpu 30 \
--chunk_size 100 \
--regions chrom2k 2000 \
--regions chrom5k 5000 \
--quantifiers chrom2k count CGN,CHN \
--quantifiers chrom5k count CGN,CHN \
--quantifiers chrom2k hypo-score CGN cutoff=0.9 \
--quantifiers chrom2k hypo-score CHN cutoff=0.9 \
--quantifiers chrom5k hypo-score CGN cutoff=0.9 \
--quantifiers chrom5k hypo-score CHN cutoff=0.9 


allcools generate-dataset  \
--allc_table AllcPaths.tsv \
--output_path GSC_epitherapy_CTCF_new1.mcds \
--chrom_size_path /scratch/devtools/mmoore/genomes/snm3c/hg38/chrom_sizes.txt \
--obs_dim cell  \
--cpu 15 \
--chunk_size 50 \
--regions CTCF /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/5_snm3c/all_CTCF_anno_sorted.bed \
--regions CTCF_500bp /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/5_snm3c/all_CTCF_500bp_anno_sorted.bed \
--regions CTCF_1k /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/5_snm3c/all_CTCF_1k_anno_sorted.bed \
--regions CTCF_scaled /scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/5_snm3c/all_CTCF_scaled_anno_sorted.bed \
--quantifiers CTCF count CGN,CHN \
--quantifiers CTCF_500bp count CGN,CHN \
--quantifiers CTCF_1k count CGN,CHN \
--quantifiers CTCF_scaled count CGN,CHN \
--quantifiers CTCF hypo-score CGN cutoff=0.9 \
--quantifiers CTCF hypo-score CHN cutoff=0.9 \
--quantifiers CTCF_500bp hypo-score CGN cutoff=0.9 \
--quantifiers CTCF_500bp hypo-score CHN cutoff=0.9 \
--quantifiers CTCF_1k hypo-score CGN cutoff=0.9 \
--quantifiers CTCF_1k hypo-score CHN cutoff=0.9 \
--quantifiers CTCF_scaled hypo-score CGN cutoff=0.9 \
--quantifiers CTCF_scaled hypo-score CHN cutoff=0.9 


hicluster prepare-impute --cell_table contact_table_rmbkl.tsv --batch_size 7467 --pad 1 --cpu_per_job 30 --chr1 1 --pos1 2 --chr2 5 --pos2 6 --output_dir impute/50K/ --chrom_size_path /scratch/devtools/mmoore/genomes/snm3c/hg38/hg38.main.chrom.sort.sizes --output_dist 500000000 --window_size 500000000 --step_size 500000000 --resolution 50000


zcat ../schicluster/CellMetadata.PassQC.Quintile.csv.gz | grep -n "B36_0" | cut -d: -f1 | awk '{print $1 - 2}' > B36_0_cells.txt
zcat ../schicluster/CellMetadata.PassQC.Quintile.csv.gz | grep -n "B36_1" | cut -d: -f1 | awk '{print $1 - 2}' > B36_1_cells.txt
zcat ../schicluster/CellMetadata.PassQC.Quintile.csv.gz | grep -n "B36_2" | cut -d: -f1 | awk '{print $1 - 2}' > B36_2_cells.txt
zcat ../schicluster/CellMetadata.PassQC.Quintile.csv.gz | grep -n "B36_3" | cut -d: -f1 | awk '{print $1 - 2}' > B36_3_cells.txt
zcat ../schicluster/CellMetadata.PassQC.Quintile.csv.gz | grep -n "B36_4" | cut -d: -f1 | awk '{print $1 - 2}' > B36_4_cells.txt
zcat ../schicluster/CellMetadata.PassQC.Quintile.csv.gz | grep -n "B36_5" | cut -d: -f1 | awk '{print $1 - 2}' > B36_5_cells.txt
zcat ../schicluster/CellMetadata.PassQC.Quintile.csv.gz | grep -n "B36_6" | cut -d: -f1 | awk '{print $1 - 2}' > B36_6_cells.txt

python higashi/Merge2Cool.py -c 10kb/config_GBM.JSON -o B36_6 -l B36_6_cells.txt -t selected -n
python higashi/Merge2Cool.py -c 10kb/config_GBM.JSON -o B36_5 -l B36_5_cells.txt -t selected -n
python higashi/Merge2Cool.py -c 10kb/config_GBM.JSON -o B36_4 -l B36_4_cells.txt -t selected -n
python higashi/Merge2Cool.py -c 10kb/config_GBM.JSON -o B36_3 -l B36_3_cells.txt -t selected -n
python higashi/Merge2Cool.py -c 10kb/config_GBM.JSON -o B36_2 -l B36_2_cells.txt -t selected -n
python higashi/Merge2Cool.py -c 10kb/config_GBM.JSON -o B36_1 -l B36_1_cells.txt -t selected -n
python higashi/Merge2Cool.py -c 10kb/config_GBM.JSON -o B36_0 -l B36_0_cells.txt -t selected -n

cd ../B36_0/
sbatch ~/tricks/Merge2Cool.sh B36_0_cells.txt 
cd ../B36_1/
sbatch ~/tricks/Merge2Cool.sh B36_1_cells.txt 
cd ../B36_2/
sbatch ~/tricks/Merge2Cool.sh B36_2_cells.txt 
cd ../B36_3/
sbatch ~/tricks/Merge2Cool.sh B36_3_cells.txt 
cd ../B36_4/
sbatch ~/tricks/Merge2Cool.sh B36_4_cells.txt 
cd ../B36_5/
sbatch ~/tricks/Merge2Cool.sh B36_5_cells.txt 


