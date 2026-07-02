#!/bin/bash
#!/usr/bin/awk
# programmer: Holden Liang
# Run it under the folder for domains
###################################################
# This is a adapted version of Vanilla pipeline 
# https://github.com/ren-lab/hic-pipeline
# Adapted by Yonghao Liang(Holden)
###################################################

#########################################################################################################################################################
## INPUT
# Arguments to bash script to decide on various parameters
MAXJOBS=12 #Cores per CPU

HIC_PATH="/scratch/yliang/Epitherapy_3D/data/hi-C/BRep" # path to hic files
HIC_FILE_EXTENSION=".hic"
DUMP_MATRIX="yes"
MATRIX_PATH=""
normalization="KR"

BIN_SIZE=25000 # resolution
BIN_SIZE_n=25k
WINDOW_SIZE_DI=80 # window length(bin_number) for DI score(usually 2,000,000bp)
WINDOW_SIZE_INS=20 # window length(bin_number) for insulation score(usually 500,000bp)

#########################################################################################################################################################

PIPE_PATH="/bar/yliang/softwares/domain_call_Bing_software"
JUICER_PATH="/bar/yliang/softwares/juicer/scripts/common"
REF_DIR="/bar/yliang/genomes"
FAI_PATH="/bar/yliang/genomes/public/hg38/hg38.fa.fai"

while [[ $# > 1 ]]
do
key="$1"
case $key in
    -j|--maxjobs)
    MAXJOBS="$2"
    shift # past argument
    ;;
    -h|--hicpath)
    HIC_PATH="$2"
    shift # past argument
    ;;
    -dm|--dumpmatrix)
    DUMP_MATRIX="$2"
    shift # past argument
    ;;
    -m|--matrixpath)
    MATRIX_PATH="$2"
    shift # past argument
    ;;
    -b|--binsize)
    BIN_SIZE="$2"
    shift # past argument
    ;;
    -bn|--binsizen)
    BIN_SIZE_n="$2"
    shift # past argument
    ;;
    --normalization)
    normalization="$2"
    shift # pass argument
    ;;
    --window_size_di)
    WINDOW_SIZE_DI="$2"
    shift # past argument
    ;;
    --window_size_ins)
    WINDOW_SIZE_INS="$2"
    shift # past argument
    ;;
    --sampletype)
    sampletype="$2"
    shift # pass argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

echo "parallel job limit: ${MAXJOBS}"
echo "hic file location: ${HIC_PATH}"
echo "resolution/bin size: ${BIN_SIZE}"
echo "Window size when computing DI: $((BIN_SIZE*WINDOW_SIZE_DI))"
echo "Window size when computing insulation score: $((BIN_SIZE*WINDOW_SIZE_INS))"

module load java matlab python3
module load R/3.6.1

echo "python3 version:" >> softwares_version.txt
python3 --version >> softwares_version.txt
echo "java version:" >> softwares_version.txt
java -version 2>> softwares_version.txt
echo "R version:" >> softwares_version.txt
R --version >> softwares_version.txt
echo "juicer version:" >> softwares_version.txt
java -jar /bar/yliang/softwares/juicer/scripts/common/juicer_tools.jar >> softwares_version.txt

echo "Resolution: $BIN_SIZE" > WindowsizeInfo.txt
echo "Window size when computing DI: $((BIN_SIZE*WINDOW_SIZE_DI))" >> WindowsizeInfo.txt
echo "Window size when computing insulation score: $((BIN_SIZE*WINDOW_SIZE_INS))" >> WindowsizeInfo.txt

mkdir RawResult
cd RawResult

## get DI score from hic files
find ${HIC_PATH} -maxdepth 1 -name "*${HIC_FILE_EXTENSION}" | while read file ; do xbase=$(basename $file ${HIC_FILE_EXTENSION}) ; echo "bash ${PIPE_PATH}/1_hic2di.sh --basename $xbase --hic_path $HIC_PATH --bin_size $BIN_SIZE --bin_size_n $BIN_SIZE_n --REF_DIR $REF_DIR --window_size_di $WINDOW_SIZE_DI --window_size_ins $WINDOW_SIZE_INS --pipe_path $PIPE_PATH --fai $FAI_PATH --dump_matrix $DUMP_MATRIX --matrix_path $MATRIX_PATH --juicer_path $JUICER_PATH --sampletype $sampletype --normalization $normalization" >> 1_hic2di_commands.txt ; done ;

parallel_GNU -j $MAXJOBS < 1_hic2di_commands.txt

## from DI to TAD
find ${HIC_PATH} -maxdepth 1 -name "*${HIC_FILE_EXTENSION}" | while read file ; do xbase=$(basename $file ${HIC_FILE_EXTENSION}) ; echo "bash ${PIPE_PATH}/2_di2tad.sh --basename $xbase --hic_path $HIC_PATH --bin_size $BIN_SIZE --bin_size_n $BIN_SIZE_n --REF_DIR $REF_DIR --window_size_di $WINDOW_SIZE_DI --window_size_ins $WINDOW_SIZE_INS --pipe_path $PIPE_PATH --fai $FAI_PATH" >> 2_di2tad_commands.txt ; done ;

parallel_GNU -j $MAXJOBS < 2_di2tad_commands.txt

## remove Y chromosome and only keep autosome + chrX
cd ../
find ./RawResult -maxdepth 1 -name "*.norm.combine.DI" | while read file ; do xbase=$(basename $file ".norm.combine.DI") ; echo "grep -P \"chrY\\t\" -v $file > $xbase.noY.DI" >> 3_CleanUp.txt ; done
find ./RawResult -maxdepth 1 -name "*k.TAD" | while read file ; do xbase=$(basename $file ".TAD") ; echo "grep -P \"chrY\\t\" -v $file > $xbase.noY.TAD" >> 3_CleanUp.txt ; done
parallel_GNU -j $MAXJOBS < 3_CleanUp.txt




