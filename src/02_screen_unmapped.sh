#!/bin/bash
#BSUB -n 8 -R "rusage[mem=8]" -R "span[hosts=1]" -W 359

## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -E -e -x -u -o pipefail

## trap read debug

echo $LSB_JOBINDEX
echo $LSB_MAX_NUM_PROCESSORS
## [[ ! -z "$LSB_MAX_NUM_PROCESSORS" ]] || exit 101

R1=( $( ls ${1}*dd.bam ) )

EXTENSION=.bam

MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//')
echo $MID

[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/

[ -d unmapped_out/${MID} ] || mkdir -p unmapped_out/${MID}
fastq_screen --aligner 'bwa' --threads $LSB_MAX_NUM_PROCESSORS --top 1000000,1000 --conf ~/conf/fastq_screen.conf --outdir unmapped_out/${MID} <( samtools fastq -f4 ${R1[$LSB_JOBINDEX]} ) --force  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastq_screen.log
if [ $? -eq 0 ] ; then echo "fastqc_screen succesfull" ; else exit 13 ; fi
