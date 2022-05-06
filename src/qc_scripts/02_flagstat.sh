#!/bin/bash
#BSUB -n 4 -R "rusage[mem=2]" -R "span[hosts=1]" -W 12:00

## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -E -e -x -u -o pipefail
echo $LSB_JOBINDEX
echo $LSB_MAX_NUM_PROCESSORS

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

OUT=$3
[[ ! -z "$OUT" ]] || OUT=bowtie_out/

ALIGNER=$(echo $OUT | sed -e 's/_.*//' -e 's:/::g')

## get BAM
BAM=$(find $OUT  -name "${MID}.*md.bam")


[ -d flagstats ] || mkdir flagstats

samtools flagstat -@ $LSB_MAX_NUM_PROCESSORS ${BAM} > flagstats/${MID}.${ALIGNER}.flagstat  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_flagstat.log
