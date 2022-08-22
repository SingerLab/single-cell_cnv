#!/bin/bash
#BSUB -n 1 -M 12 -W 48:00
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This script to runs a de-multiplexing tools to generate paired-end"
    echo " (R1 and R2) reads from each cell."
    echo ""
    echo "Usage:"    
    echo "This script expects the FASTQ with all cells as the first argument, "
    echo " the FASTQ extension as the second e.g. _001.fastq.gz, .fq.gz;"
    echo " the output directory name as the third argument"
    echo ""
    echo "Example:"
    echo "bsub -n 1 -M 12 -W 89 ./src/bdplex/barcode.wsplit.sh path/to/Sample_R1_001.fastq.gz _R1_001.fastq.gz wsplit/"
    echo "bsub -J 'wsplit[1-1500]' -n 1 -M 12 -W 356 ./src/bdplex/barcode.wsplit.sh path/to/Sample_R1_001.fastq.gz _R1_001.fastq.gz wsplit/"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 1 -M 12 -W 89 LSB_JOBINDEX=9 ./src/bdplex/barcode.wsplit.sh path/to/Sample_R1_001.fastq.gz _R1_001.fastq.gz wsplit/"
    echo ""
    exit 1;
}
#</usage>
set -e -x -o pipefail -u

R1=$1
R2=$( echo ${R1} | sed -e 's/_R1/_R2/' )

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fastq.gz

X1=$( basename $R1 $EXTENSION )

BWGA192_IDX=( $( cut -f 1 src/bdeplex/index.192.txt ) )
BWGA192_SEQ=( $( cut -f 2 src/bdeplex/index.192.txt ) )

echo "Searching for" ${BWGA192_IDX[$LSB_JOBINDEX]} ${BWGA192_SEQ[$LSB_JOBINDEX]}

OUT=$3
[[ ! -z "$OUT" ]] || OUT=wsplit/
[ -d $OUT ] || mkdir -p $OUT

## if nreads file exists: exit
[[ ! -f ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.nreads.txt ]] || exit 1

## use of wannaAln
wannaAln -a ${R1} -b ${R2} -x $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz \
	 -y $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz \
	 -q ${BWGA192_SEQ[$LSB_JOBINDEX]} -m 0
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} barcode split complete" > $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi

echo ${X1}.${BWGA192_IDX[$LSB_JOBINDEX]} \
     $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz | awk 'NR%4==1' | wc -l )  \
     $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz | awk 'NR%4==1' | wc -l ) | \
    tr ' ' "\t" > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.nreads.txt
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} read count complete" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi
