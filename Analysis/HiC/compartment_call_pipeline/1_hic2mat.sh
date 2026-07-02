#!/bin/bash
#!/usr/bin/awk
#####################################################################
# This is a script for A/B compartments calling using Matlab.
# Adapted from Dixon, Jesse R., et al. 2015, by Xianglin Zhang.
# Input is .hic files, output is dumped matrix
#####################################################################
while [[ $# > 1 ]]
do
	key="$1"
	case $key in
		--sample)
		sample="$2"
		shift # pass argument
		;;
		--basename)
		basename="$2"
		shift # pass argument
		;;
		--matrix_path)
		matrix_path="$2"
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
		--normalization)
		normalization="$2"
		shift # pass argument
		;;
		--REF_DIR)
		REF_DIR="$2"
		shift # pass argument
		;;
		--juicer_path)
		juicer_path="$2"
		shift # pass argument
    	;;
    	--pipe_path)
		pipe_path="$2"
		shift
		;;
    	*)
    	        # unknown option
    	;;
	esac
	shift # past argument or value
done

ml java

for chr in `awk '{print $1}' $REF_DIR/private/hg38_23_chromosomes.simplified.size|sort -n`; do
	chr_size=$(awk -v chr=$chr '{{if($1==chr)print $3}}' $REF_DIR/private/hg38_23_chromosomes.simplified.size)
	java -jar $juicer_path/juicer_tools.jar dump observed ${normalization} ${sample} $chr $chr BP ${bin_size} |\
		awk '{if (NR!=1) {print $0}}' |\
			$pipe_path/col2mat.awk -v chr=$chr -v bin_size=${bin_size} -v chr_size=$chr_size |\
				awk '{if (NR!=1) {print $0}}'\
					> ${matrix_path}/${basename}.${bin_size_n}.chr$chr.${normalization}.observed.matrix # in col2mat.awk, NaN needs to be kept
	if [ $chr == "X" ]
	then
		mv ${matrix_path}/${basename}.${bin_size_n}.chrX.${normalization}.observed.matrix ${matrix_path}/${basename}.${bin_size_n}.chr23.${normalization}.observed.matrix
	fi
done
