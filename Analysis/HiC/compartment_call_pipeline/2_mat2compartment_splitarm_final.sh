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
		--matrix_path)
		matrix_path="$2"
		shift # pass argument
		;;
		--output_path)
		output_path="$2"
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

export basename=${basename}
export prefix=$matrix_path/${basename}.${bin_size_n}.chr
export suffix=.${normalization}.observed.matrix
export storage=$output_path
export pipe_path=$pipe_path
export bin_size=$bin_size
#matlab -nodisplay -nosplash -nodesktop < $pipe_path/compartment_call_nan_for_no_signal_bins.m > $output_path/${basename}.$bin_size_n.bing2015.compartment_call.log 2> $output_path/${basename}.$bin_size_n.bing2015.compartment_call.err
matlab -nodisplay -nosplash -nodesktop < $pipe_path/compartment_call_nan_for_no_signal_bins.splitarm_final.m > $output_path/${basename}.$bin_size_n.bing2015.compartment_call.log 2> $output_path/${basename}.$bin_size_n.bing2015.compartment_call.err
awk '{print "chr"$1,$2,$3,$4}' $output_path/${basename}.AB_compartment_4_0_single_arm.txt | awk -v OFS="\t" '$1=$1' | awk -v OFS="\t" '{if($1=="chr23") print("chrX",$2,$3,$4); else print $0;}'  > $output_path/final.${basename}.bing2015.eigen.${normalization}.${bin_size_n}.bedgraph


