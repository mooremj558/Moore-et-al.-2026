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
		--juicer_path)
		juicer_path="$2"
		shift # pass argument
		;;
		--normalization)
		normalization="$2"
		shift # pass argument
		;;
		--matrix_path)
		matrix_path="$2"
		shift # pass argument
		;;
		--dump_matrix)
		dump_matrix="$2"
		shift # pass argument
		;;
		--sampletype)
		sampletype="$2"
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
for chr in `awk '{print $2}' $REF_DIR/private/hg38_23_chromosomes.simplified.size|sort -V`; do
	chr_size=$(awk -v chr=$chr '{{if($2==chr)print $3}}' $REF_DIR/private/hg38_23_chromosomes.simplified.size)
	### .hic to matrix using straw
	#echo "	$pipe_path/scripts/straw ${normalization} ${hic_path}/${basename}.hic $chr $chr BP $bin_size 
	#	$pipe_path/scripts/col2mat.awk -v chr=$chr -v bin_size=$bin_size -v chr_size=$chr_size 
	#		awk '{if (NR!=1) {print $0}}'" > ${basename}.${bin_size_n}.matrix_calls.log
	if [ $dump_matrix == "yes" ]; then
		echo "generating matrix"
		#java -jar $juicer_path/juicer_tools.jar dump observed ${normalization} ${hic_path}/${basename}.hic $chr $chr BP ${bin_size} |\
		#	awk 'NR!=1' |\
		#		$pipe_path/scripts/col2mat.awk -v chr=$chr -v bin_size=${bin_size} -v chr_size=$chr_size \
		#			> ${matrix_path}/${basename}.${bin_size_n}.chr$chr.${normalization}.observed.matrix
		if [[ $sampletype == "our" ]]; then
			starw_chr=${chr#"chr"}
			# /bar/yliang/softwares/domain_call_Bing_software/scripts/straw ${normalization} /scratch/yliang/Epitherapy_3D/data/hi-C/Mega/B36_DP.hic 5 5 BP 10000
			# java -jar /bar/yliang/softwares/juicer/scripts/common/juicer_tools.jar dump observed ${normalization} /scratch/yliang/Epitherapy_3D/data/hi-C/Mega/B36_DP.hic 5 5 BP 10000|awk 'NR!=1' |\
			$pipe_path/scripts/straw ${normalization} ${hic_path}/${basename}.hic $starw_chr $starw_chr BP $bin_size |\
				$pipe_path/scripts/col2mat.awk -v chr=$chr -v bin_size=$bin_size -v chr_size=$chr_size |\
					awk '{if (NR!=1) {print $0}}'\
						> ${matrix_path}/${basename}.${bin_size_n}.$chr.${normalization}.observed.matrix 2> ${basename}.${bin_size_n}.matrix_calls.err # straw will have nan in the matrix, need to convert them into 0. If use juicer dump to get matrix, there will be NaN in the matrix, this awk converts both to 0
		else
			$pipe_path/scripts/straw ${normalization} ${hic_path}/${basename}.hic $chr $chr BP $bin_size |\
				$pipe_path/scripts/col2mat.awk -v chr=$chr -v bin_size=$bin_size -v chr_size=$chr_size |\
					awk '{if (NR!=1) {print $0}}'\
						> ${matrix_path}/${basename}.${bin_size_n}.$chr.${normalization}.observed.matrix 2> ${basename}.${bin_size_n}.matrix_calls.err
		fi
		echo "${basename}.${bin_size_n}.$chr.${normalization}.observed.matrix is done"
	else
		echo "seems like you have already generated matrix"
	fi
	### matrix to directionality index
	echo "generating DI"
	Rscript $pipe_path/scripts/asc2di.R ${matrix_path}/${basename}.${bin_size_n}.$chr.${normalization}.observed.matrix $chr $bin_size $window_size_di ${basename}.${bin_size_n}.$chr.norm.DI
	echo "${basename}.${bin_size_n}.$chr.norm.DI is done"
	### matrix to insulation score
	Rscript $pipe_path/scripts/mat2insulation.R -m ${matrix_path}/${basename}.${bin_size_n}.$chr.${normalization}.observed.matrix -b $bin_size -w $((bin_size*window_size_ins)) -c $chr -o ${basename}.${bin_size_n}.$chr.insulation.bed 2> ${basename}.${bin_size_n}.insulation_calls.err
	cat ${basename}.${bin_size_n}.${chr}.norm.DI >> ${basename}.${bin_size_n}.norm.combine.DI
	cat ${basename}.${bin_size_n}.${chr}.insulation.bed >> ${basename}.${bin_size_n}.combine.insulation.bed
	#rm ${basename}.${bin_size_n}.chr${chr}.norm.asc
	rm ${basename}.${bin_size_n}.${chr}.norm.DI
	rm ${basename}.${bin_size_n}.${chr}.insulation.bed
done

sed 's/chr//g' ${basename}.${bin_size_n}.norm.combine.DI | sed 's/X/23/g' | grep -v Y > ${basename}.${bin_size_n}.norm.combine.forHMM.DI 

#awk -v OFS="\t" '{if($1=="X"){print 23,$2,$3,$4} else if($1=="Y"){print 24,$2,$3,$4} else{print($1,$2,$3,$4)}}' ${basename}.${bin_size_n}.norm.combine.DI > ${basename}.${bin_size_n}.norm.combine.modified.DI

#awk -v OFS="\t" '{print("chr"$1,$2,$3,$4)}' ${basename}.${bin_size_n}.norm.combine.DI > ${basename}.${bin_size_n}.norm.combine.for_browser.DI

