#!/bin/bash
#BSUB -n 1 -R "rusage[mem=6]" -R "span[hosts=1]" -W 12:00

## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -E -e -x -u -o pipefail

echo $LSB_JOBINDEX
echo $LSB_MAX_NUM_PROCESSORS

R1=( $( ls ${1}*fq.gz ) )

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fq.gz

MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//')
echo $MID

[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/

[ -d fastp/ ] || mkdir -p fastp/

fastp -i ${R1[$LSB_JOBINDEX]} -p -h fastp/${MID}.html -j fastp/${MID}.json -R ${MID} --umi --umi_loc read1 --umi_len 6  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastp.log
if [ $? -eq 0 ] ; then echo "fastp succesfull" ; touch fastp/${MID}.fastp.ok ; else exit 11 ; fi

