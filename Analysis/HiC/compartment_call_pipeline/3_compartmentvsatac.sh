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
		--basename)
		basename="$2"
		shift # pass argument
		;;
		--bin_size_n)
		bin_size_n="$2"
		shift # pass argument
		;;
		--peak_path)
		peak_path="$2"
		shift
		;;
    	--pipe_path)
		pipe_path="$2"
		shift
		;;
		--normalization)
		normalization="$2"
		shift
		;;
		--REF_DIR)
		REF_DIR="$2"
		shift # pass argument
		;;
    	*)
    	        # unknown option
    	;;
	esac
	shift # past argument or value
done

sample="final.${basename}.bing2015.eigen.${normalization}.${bin_size_n}.bedgraph"
peaksignal="${basename}.ATAC_signal.compartment_call_bing_2015.${bin_size_n}.peaksignal"
printf "\n$basename\n" >> 3_compartmentvsatac.summary.txt
export input=$sample
export peaknumber=$peak_path/$peaksignal
export output=3_compartmentvsatac.summary.txt
export split=$REF_DIR/private/hg38_chrom_split.txt 
python3 $pipe_path/3_compartmentvsatac_${bin_size_n}_resolution.py 2> 3_compartmentvsatac.summary.txt.err




