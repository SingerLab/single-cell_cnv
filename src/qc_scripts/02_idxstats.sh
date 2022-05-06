#!/bin/bash
#BSUB -n 4 -R "rusage[mem=2]" -R "span[hosts=1]" -W 12:00

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

OUT=$3
[[ ! -z "$OUT" ]] || OUT=bowtie_out/

ALIGNER=$(echo $OUT | sed -e 's/_.*//' -e 's:/::g')

## consider using something like this to assign aligner
## VAR='GNU/Linux is an operating system'
## if [[ $VAR =~ .*Linux.* ]]; then
##   echo "It's there."
## fi

## get BAM
BAM=$(find $OUT  -name "${MID}.*md.bam")


[ -d idxstats ] || mkdir idxstats

samtools idxstats -@ $LSB_MAX_NUM_PROCESSORS ${BAM} > idxstats/${MID}.${ALIGNER}.idxstats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log


