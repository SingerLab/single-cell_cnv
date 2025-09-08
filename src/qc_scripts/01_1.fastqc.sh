#!/bin/bash
#BSUB -n 8 -R "rusage[mem=1]" -R "span[hosts=1]" -W 12:00

## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -E -e -x -u -o pipefail

## trap read debug

MAX_MEM_GB=32G

echo $LSB_JOBINDEX
echo $LSB_MAX_NUM_PROCESSORS
## [[ ! -z "$LSB_MAX_NUM_PROCESSORS" ]] || exit 101

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fq.gz

R1=( $( ls ${1}*${EXTENSION} ) )


MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//')
echo $MID

[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/

[ -d fastqc/ ] || mkdir -p fastqc/

fastqc -t $LSB_MAX_NUM_PROCESSORS -a ~/genomes/adapters/bwga/adapter_list.txt -o fastqc/ ${R1[$LSB_JOBINDEX]}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastqc.log

