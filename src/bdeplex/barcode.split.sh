#!/bin/bash
#BSUB -n 1 -M 2 -W 359
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This script demultiplexes an R1 FASTQ file based on a preset array of"
    echo " 192 barcodes. One --user specified-- barcode is demupltiplexed. Barcodes"
    echo " are specified via the `LSB_JOBINDEX` environment variable.  Script also"
    echo " counts the fastq file and provides an `nreads.txt`"
    echo ""
    echo "Usage:"    
    echo "This script expects a directory name with the FASTQ files as the first"
    echo " argument; the FASTQ extension as the second e.g. _001.fastq.gz, .fq.gz;"
    echo " the output directory name as the third argument."
    echo ""
    echo "Example:"
    echo "bsub -n 1 -M 2 -W 89 ./src/bdeplex/barcode.split.sh path/to/file.fastq.gz .fastq.gz bsplit/"
    echo ""
    echo "bsub -J 'bm[1-191]' -n 1 -M 2 -W 89 ./src/bdeplex/barcode.split.sh path/to/file.fastq.gz .fastq.gz bsplit/"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 8 -M 4 -W 89 LSB_JOBINDEX=9  ./src/bdeplex/barcode.split.sh path/to/file.fastq.gz .fastq.gz bsplit/"
    echo ""
    exit 1;
}
#</usage>

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

