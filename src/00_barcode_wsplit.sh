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
    echo "bsub -n 1 -M 2 -W 89 ./src/bdeplex/barcode.split.sh path/to/file_R1.fastq.gz path/to/file_R2.fastq.gz _R1.fastq.gz wsplit/"
    echo ""
    echo "bsub -J 'bm[1-191]' -n 1 -M 2 -W 89 ./src/bdeplex/barcode.split.sh path/to/file_R1.fastq.gz path/to/file_R2.fastq.gz _R1.fastq.gz wsplit/"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 8 -M 4 -W 89 LSB_JOBINDEX=9  ./src/bdeplex/barcode.split.sh path/to/file_R1.fastq.gz path/to/file_R2.fastq.gz _R1.fastq.gz wsplit/"
    echo ""
    exit 1;
}
#</usage>

set -e -x -o pipefail -u

R1=$1
R2=$2

EXTENSION=$3
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fastq.gz

X1=$( basename $R1 $EXTENSION )

BWGA192_IDX=( $( cut -f 1 src/bdeplex/barcode.192.txt ) )
BWGA192_SEQ=( $( cut -f 2 src/bdeplex/barcode.192.txt ) )

echo "searching for" ${BWGA192_IDX[$LSB_JOBINDEX]} ${BWGA192_SEQ[$LSB_JOBINDEX]}

OUT=$4
[[ ! -z "$OUT" ]] || OUT=wsplit/
[ -d $OUT ] || mkdir $OUT

module load gcc/11.2.0

BARCODE_SPLITTER=/work/singer/opt/src/barcode_split/barcode_split
CHECK_FASTQ_PAIRS=/work/singer/opt/src/check_fastq_pairs/check_fastq_pairs

## OUTPUT FILE STEM
MID=${OUT}/${X1}_${BWGA192_IDX[$LSB_JOBINDEX]}_${BWGA192_IDX[$LSB_JOBINDEX]}

## checkpoint
OK_FILE=${MID}.ok

# Function to update checkpoint in OK file
update_checkpoint() {
    echo "## $(date) $1 complete" >> "$OK_FILE"
}

# Function to check if step has been completed
check_checkpoint() {
    grep -q "$1 complete" "$OK_FILE" 2> /dev/null
}

# Pipeline start
if [ ! -f $OK_FILE ] ; then
    echo "## $(date) wsplit start" > "$OK_FILE"
fi

if ! check_checkpoint "barcode split"; then
    ${BARCODE_SPLITTER} -p ${BWGA192_SEQ[$LSB_JOBINDEX]} -o ${MID} ${R1} ${R2}
    if [ $? -eq 0 ] ; then
	update_checkpoint "barcode split"
    else
	exit 2
    fi
fi

## check barcode specific R1/R2
if ! check_checkpoint "barcode fq_check"; then
${CHECK_FASTQ_PAIRS} ${MID}_R{1,2}.fastq.gz > \
		     ${MID}.fq_check.txt
    if [ $? -eq 0 ] ; then
	update_checkpoint "barcode fq_check"
    else
	exit 2
    fi
fi


#__EOF__
