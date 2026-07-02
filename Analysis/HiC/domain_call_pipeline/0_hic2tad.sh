#!/bin/bash
#!/usr/bin/awk
###################################################
# This is a adapted version of Vanilla pipeline 
# Adapted from Dixon, Jesse R., et al. 2015
# from Bing Ren's lab to process Hi-C data to .TAD
# https://github.com/ren-lab/hic-pipeline
# Adapted by Yonghao Liang(Holden)
###################################################
while [[ $# > 1 ]]
do
	key="$1"
	case $key in
		--basename)
		basename="$2"
		shift # pass argument
		;;
		--bin_size)
		bin_size="$2"
		shift # pass argument
		;;
		--bin_size_n)
		bin_size_n="$2"
		shift # pass argument
		;;
		--window_size_di)
		window_size_di="$2"
		shift # pass argument
		;;
		--window_size_ins)
		window_size_ins="$2"
		shift # pass argument
		;;
		--hic_path)
		hic_path="$2"
		shift # pass argument
		;;
		--pipe_path)
		pipe_path="$2"
		shift # pass argument
		;;
		--REF_DIR)
		REF_DIR="$2"
		shift
		;;
		--fai)
		fai="$2"
		shift
		;;
    	*)
    	        # unknown option
    	;;
	esac
	shift # past argument or value
done

touch ${basename}.${bin_size_n}.combine.insulation.bed
touch ${basename}.${bin_size_n}.norm.combine.DI
for chr in `awk '{print $1}' $REF_DIR/private/hg38_23_chromosomes.simplified.size|sort -n`; do
	chr_size=$(awk -v chr=$chr '{{if($1==chr)print $3}}' $REF_DIR/private/hg38_23_chromosomes.simplified.size)
	### .hic to matrix using straw
	$pipe_path/scripts/straw KR $hic_path/${basename}.hic $chr $chr BP $bin_size |\
		$pipe_path/scripts/col2mat.awk -v chr=$chr -v bin_size=$bin_size -v chr_size=$chr_size \
			> ${basename}.${bin_size_n}.chr${chr}.norm.asc
	echo "${basename}.${bin_size_n}.chr${chr}.norm.asc is done"
	### matrix to directionality index
	Rscript $pipe_path/scripts/asc2di.R ${basename}.${bin_size_n}.chr$chr.norm.asc $chr $bin_size $window_size_di ${basename}.${bin_size_n}.chr$chr.norm.DI
	echo "${basename}.${bin_size_n}.chr$chr.norm.DI is done"
	### matrix to insulation score
	Rscript $pipe_path/scripts/mat2insulation.R -m ${basename}.${bin_size_n}.chr$chr.norm.asc -b $bin_size -w $((bin_size*window_size_ins)) -c $chr -o ${basename}.${bin_size_n}.chr$chr.insulation.bed
	cat ${basename}.${bin_size_n}.chr${chr}.norm.DI >> ${basename}.${bin_size_n}.norm.combine.DI
	cat ${basename}.${bin_size_n}.chr${chr}.insulation.bed >> ${basename}.${bin_size_n}.combine.insulation.bed
	rm ${basename}.${bin_size_n}.chr${chr}.norm.asc
	rm ${basename}.${bin_size_n}.chr${chr}.norm.DI
	rm ${basename}.${bin_size_n}.chr${chr}.insulation.bed
done

awk '{if($1=="X"){print 23,$2,$3,$4} else if($1=="Y"){print 24,$2,$3,$4} else{print$0}}' ${basename}.${bin_size_n}.norm.combine.DI > ${basename}.${bin_size_n}.norm.combine.modified.DI
### DI to TAD using domaincall software
echo "begin HMM_calls"
export HMM_input=${basename}.${bin_size_n}.norm.combine.modified.DI
export HMM_output=${basename}.${bin_size_n}.norm.combine.hmm.DI
export HMM_PATH=$pipe_path
module load matlab
matlab -nodisplay -nosplash -nodesktop < $pipe_path/HMM_calls.m > $basename.${bin_size_n}.HMM_calls.log 2> ${basename}.${bin_size_n}.HMM_calls.err
echo "finish HMM_calls"
perl $pipe_path/scripts/file_ends_cleaner.pl ${basename}.${bin_size_n}.norm.combine.hmm.DI ${basename}.${bin_size_n}.norm.combine.modified.DI | perl $pipe_path/scripts/converter_7col.pl | sed 's/chrchr/chr/' > $basename.${bin_size_n}.norm.combine.hmm.7col
echo "begin to generate .TAD file"
min=2; prob=0.99
for chr in `awk '{print $1}' ${basename}.${bin_size_n}.norm.combine.hmm.7col|sort -u`; do
	#perl $pipe_path/scripts/hmm_probablity_correcter.pl <(awk -v chr=$chr '{if($1==chr) print $0}' ${basename}.${bin_size_n}.norm.combine.hmm.7col) $min $prob $bin_size | perl $TAD_PERL/hmm-state_caller.pl $fai $chr | perl $TAD_PERL/hmm-state_domains.pl
	perl $pipe_path/scripts/hmm_probablity_correcter.pl <(awk '/^'$chr'\t/' ${basename}.${bin_size_n}.norm.combine.hmm.7col) $min $prob $bin_size | perl $pipe_path/scripts/hmm-state_caller.pl $fai $chr | perl $pipe_path/scripts/hmm-state_domains.pl
done > final.${basename}.${bin_size_n}.TAD
echo "all is done"
	
