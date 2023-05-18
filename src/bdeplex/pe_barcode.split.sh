#!/bin/bash
#BSUB -n 1 -M 16 -W 24:00
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "Identifying corcodant barcodes in R1 and R2 files"
    echo "This script requires to run the barcode.wsplit.sh first"
    exit 1;
}
#</usage>
set -e -x -o pipefail -u

R1=$1
R2=$( echo ${R1} | sed -e 's/R1/R2/' )

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fq.gz

X1=$( basename $R1 $EXTENSION )

BC=$( echo $X1 | sed -e 's/.*b//' )
echo "Barcode is" $BC

if [ $BC -eq 1 ]
then
    BCIDX=0
elif [ $BC -ne 1 ]
then
    BCIDX=$(expr $BC - 1)
fi

echo $BCIDX

BWGA192_IDX=( $( cut -f 1 src/bdeplex/index.192.txt ) )
BWGA192_SEQ=( $( cut -f 2 src/bdeplex/index.192.txt ) )

echo "Searching for" ${BWGA192_IDX[$BCIDX]} ${BWGA192_SEQ[$BCIDX]}

OUT=$3
[[ ! -z "$OUT" ]] || OUT=wsplit2/
[ -d $OUT ] || mkdir -p $OUT/

tmp_out_file=${OUT}/tmp.${X1}.${BWGA192_IDX[$BCIDX]}.R1.R2.txt.gz

R1_OUT_FILE=${OUT}/${X1}.${BWGA192_IDX[$BCIDX]}.R1.fq.gz
R2_OUT_FILE=${OUT}/${X1}.${BWGA192_IDX[$BCIDX]}.R2.fq.gz

paste <( zcat ${R1} | paste  - - - - )  \
      <( zcat ${R2} | paste  - - - - ) | \
    awk '$2 ~ /'${BWGA192_SEQ[$BCIDX]}'[ACTGN]{83,}/ &&
         $6 ~ /'${BWGA192_SEQ[$BCIDX]}'[ACTGN]{83,}/ {print}' | \
    gzip -c > ${tmp_out_file}

zcat ${tmp_out_file} | cut -f 1-4 | tr "\t" "\n" | gzip -c > $R1_OUT_FILE
zcat ${tmp_out_file} | cut -f 5-8 | tr "\t" "\n" | gzip -c > $R2_OUT_FILE

rm ${tmp_out_file}

echo ${X1}.${BWGA192_IDX[$BCIDX]} \
     $( zcat ${R1_OUT_FILE} | awk 'NR%4==1' | wc -l )  \
     $( zcat ${R2_OUT_FILE} | awk 'NR%4==1' | wc -l ) | \
    tr ' ' "\t" > ${OUT}/${X1}.${BWGA192_IDX[$BCIDX]}.nreads.txt



#% ## old script
#% set -e -x
#% 
#% ## conda init
#% ## conda activate fastq-manip
#% 
#% R1=$1
#% R2=$2
#% 
#% EXTENSION=$3
#% [[ ! -z "$EXTENSION" ]] || EXTENSION=.fastq.gz
#% 
#% X1=$( basename $R1 $EXTENSION )
#% 
#% BWGA192_IDX=( $( cut -f 1 src/bdeplex/barcode.192.txt ) )
#% BWGA192_SEQ=( $( cut -f 2 src/bdeplex/barcode.192.txt ) )
#% 
#% echo "searching for: " ${BWGA192_IDX[$LSB_JOBINDEX]} ${BWGA192_SEQ[$LSB_JOBINDEX]}
#% 
#% OUT=$4
#% [[ ! -z "$OUT" ]] || OUT=bsplit/
#% [ -d $OUT ] || mkdir $OUT
#% 
#% ## zgrep -B 1 -A 2 -E "${BWGA192_SEQ[$LSB_JOBINDEX]}[ACTGN]{83}" ${R1} | sed -r '/^--$/d' | gzip -c > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz
#% if [ $? -eq 0 ] ; then echo "R1: ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz complete" ; else exit 1 ; fi
#% 
#% 
#% [ -d tmp ] || mkdir tmp
#% TMPDIR=tmp/
#% 
#% ## zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz | awk 'NR%4==1' > ${TMPDIR}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.id.list.txt
#% if [ $? -eq 0 ] ; then echo "getting list of read IDs" ; else exit 2 ; fi
#% 
#% ## disfunctional !! need alternative 
#% seqtk subseq $R2 ${TMPDIR}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.id.list.txt | gzip -c > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz
#% 
#% if [ $? -eq 0 ] ; then echo "R2: ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz complete" ; rm ${TMPDIR}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.id.list.txt ; else exit 3 ; fi

