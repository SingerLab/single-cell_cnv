#!/bin/bash
#BSUB -n 1 -M 2 -W 359
set -e -x -o pipefail -u

R1=$1

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fastq.gz

X1=$( basename $R1 $EXTENSION )

BWGA192_IDX=( $( cut -f 1 src/bdeplex/barcode.192.txt ) )
BWGA192_SEQ=( $( cut -f 2 src/bdeplex/barcode.192.txt ) )

echo "searching for" ${BWGA192_IDX[$LSB_JOBINDEX]} ${BWGA192_SEQ[$LSB_JOBINDEX]}

OUT=$3
[[ ! -z "$OUT" ]] || OUT=bsplit/
[ -d $OUT ] || mkdir $OUT

zgrep -B 1 -A 2 -E "${BWGA192_SEQ[$LSB_JOBINDEX]}[ACTGN]{83}" ${R1} | sed -r '/^--$/d' | gzip -c > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.fq.gz

echo ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]} $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.fq.gz | awk 'NR%4==1' | wc -l ) | tr ' ' "\t" > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.nreads.txt

