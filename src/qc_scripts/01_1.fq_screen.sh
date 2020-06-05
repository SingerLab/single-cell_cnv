#!/bin/bash
#BSUB -n 8 -R "rusage[mem=6]" -R "span[hosts=1]" -W 24:00

## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -E -e -x -u -o pipefail

## trap read debug

MAX_MEM_GB=32G

echo $LSB_JOBINDEX
echo $LSB_MAX_NUM_PROCESSORS
## [[ ! -z "$LSB_MAX_NUM_PROCESSORS" ]] || exit 101

BWA_h37=/ifs/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
GFF_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Annotation/Genes/gencode.v19.annotation.gtf
BOWTIE_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/BowtieIndex/b37.fasta

R1=( $( ls ${1}*fq.gz ) )

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fq.gz

MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//')
echo $MID

[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/

[ -d fastq_screen ] || mkdir -p fastq_screen/

fastq_screen --aligner bwa --threads $LSB_MAX_NUM_PROCESSORS --subset 0 --conf ~/conf/fastq_screen.conf --outdir fastq_screen/ ${R1[$LSB_JOBINDEX]}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastq_screen.log
if [ $? -eq 0 ] ; then echo "fastqc_screen succesfull" ; else exit 13 ; fi
