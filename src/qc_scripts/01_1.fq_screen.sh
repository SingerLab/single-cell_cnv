#!/bin/bash
#BSUB -n 8 -R "rusage[mem=6]" -R "span[hosts=1]" -W 24:00

## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -E -e -x -u -o pipefail

## trap read debug
MEM_PER_JOB=$(echo $(printf "%d\n" ${LSB_CG_MEMLIMIT} )/1024^3 | bc )
MAX_MEM_GB=$(( MEM_PER_JOB*LSB_MAX_NUM_PROCESSORS - 1 ))G

echo $LSB_JOBINDEX
echo $LSB_MAX_NUM_PROCESSORS
## [[ ! -z "$LSB_MAX_NUM_PROCESSORS" ]] || exit 101

FASTQ_DIR=$1
EXTENSION=$2
OUTPUT_DIR=$3

## get R1 data
R1=( $( ls ${FASTQ_DIR}*${EXTENSION} ) )

MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//')
echo $MID

[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/

[ -d $OUTPUT_DIR ] || mkdir -p $OUTPUT_DIR

fastq_screen --aligner bwa --threads $LSB_MAX_NUM_PROCESSORS \
	     --subset 0 \
	     --conf ~/conf/fastq_screen.conf --outdir ${OUTPUT_DIR}  ${R1[$LSB_JOBINDEX]}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastq_screen.log


	     
if [ $? -eq 0 ] ; then
    echo "fastqc_screen succesfull" ;
else
    exit 13 ;
    rm ${OUTPUT_DIR}/${R1[$LSB_JOBINDEX]}_temp_subset.fastq ;
    rm ${OUTPUT_DIR}/${MID}_screen.txt ;
fi
