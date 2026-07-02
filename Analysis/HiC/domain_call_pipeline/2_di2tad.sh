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

### DI to TAD using domaincall software
echo "begin HMM_calls"
export HMM_input=${basename}.${bin_size_n}.norm.combine.forHMM.DI
export HMM_output=${basename}.${bin_size_n}.norm.combine.hmm.DI
export HMM_PATH=$pipe_path
module load matlab
#matlab -nodisplay -nosplash -nodesktop < $pipe_path/HMM_calls.m > ${basename}.${bin_size_n}.HMM_calls.log 2> ${basename}.${bin_size_n}.HMM_calls.err
matlab -nodisplay -nosplash -nodesktop < $pipe_path/HMM_calls.new.m > ${basename}.${bin_size_n}.HMM_calls.log 2> ${basename}.${bin_size_n}.HMM_calls.err
echo "finish HMM_calls"
perl $pipe_path/scripts/file_ends_cleaner.pl ${basename}.${bin_size_n}.norm.combine.hmm.DI ${basename}.${bin_size_n}.norm.combine.forHMM.DI  | perl $pipe_path/scripts/converter_7col.pl | sed 's/chrchr/chr/' > $basename.${bin_size_n}.norm.combine.hmm.7col
echo "begin to generate .TAD file"
min=2; prob=0.99
for chr in `awk '{print $1}' ${basename}.${bin_size_n}.norm.combine.hmm.7col|sort -u`; do
	#perl $pipe_path/scripts/hmm_probablity_correcter.pl <(awk -v chr=$chr '{if($1==chr) print $0}' ${basename}.${bin_size_n}.norm.combine.hmm.7col) $min $prob $bin_size | perl $TAD_PERL/hmm-state_caller.pl $fai $chr | perl $TAD_PERL/hmm-state_domains.pl
	perl $pipe_path/scripts/hmm_probablity_correcter.pl <(awk '/^'$chr'\t/' ${basename}.${bin_size_n}.norm.combine.hmm.7col) $min $prob $bin_size | perl $pipe_path/scripts/hmm-state_caller.pl $fai $chr | perl $pipe_path/scripts/hmm-state_domains.pl
done > ${basename}.${bin_size_n}.raw.TAD

awk -v OFS="\t" '{if(NF==3) print($1,$2,$3); else if(NF==5) print($3,$4,$5)}' ${basename}.${bin_size_n}.raw.TAD | sed 's/chr23/chrX/g' > ${basename}.${bin_size_n}.TAD

paste ${basename}.${bin_size_n}.TAD ${basename}.${bin_size_n}.TAD > ${basename}.${bin_size_n}.for_juicer.TAD

echo "all is done"

