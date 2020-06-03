 #!/bin/bash
#BSUB -n 1 -M 16 -W 24:00
set -e -x

## conda init
## conda activate fastq-manip

R1=$1
R2=$2

EXTENSION=$3
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fastq.gz

X1=$( basename $R1 $EXTENSION )

BWGA192_IDX=( $( cut -f 1 src/bdeplex/barcode.192.txt ) )
BWGA192_SEQ=( $( cut -f 2 src/bdeplex/barcode.192.txt ) )

echo "searching for: " ${BWGA192_IDX[$LSB_JOBINDEX]} ${BWGA192_SEQ[$LSB_JOBINDEX]}

OUT=$4
[[ ! -z "$OUT" ]] || OUT=bsplit/
[ -d $OUT ] || mkdir $OUT

## zgrep -B 1 -A 2 -E "${BWGA192_SEQ[$LSB_JOBINDEX]}[ACTGN]{83}" ${R1} | sed -r '/^--$/d' | gzip -c > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz
if [ $? -eq 0 ] ; then echo "R1: ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz complete" ; else exit 1 ; fi


[ -d tmp ] || mkdir tmp
TMPDIR=tmp/

## zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz | awk 'NR%4==1' > ${TMPDIR}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.id.list.txt
if [ $? -eq 0 ] ; then echo "getting list of read IDs" ; else exit 2 ; fi

## disfunctional !! need alternative 
seqtk subseq $R2 ${TMPDIR}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.id.list.txt | gzip -c > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz

if [ $? -eq 0 ] ; then echo "R2: ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz complete" ; rm ${TMPDIR}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.id.list.txt ; else exit 3 ; fi

