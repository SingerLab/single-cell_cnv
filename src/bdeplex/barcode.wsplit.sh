#!/bin/bash
#BSUB -n 1 -M 6 -W 128:00
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This script to runs a de-multiplexing tools to generate paired-end"
    echo " (R1 and R2) reads from each cell."
    echo ""
    echo "Details:"
    echo "BWGA uses combinatorial indexing to identify each cell.  During the"
    echo " Klenow Fragment reaction, a BWGA universal adapter [TGTGTTGGGTGTGTTTGG]"
    echo " is inserted into each fragment.  In the subsequent DOP-PCR, a 5' tailed"
    echo " universal adapter is used to prime the reaction, and simultaenously "
    echo " introduces a cell-index to identify each cell on the plate.  The resulting"
    echo " library contains the pattern '^[ATCG]*[ATCG]{12}TGTGTTGGGTGTGTTTGG' in"
    echo " the first 50 bases of each read (R1/R2).  A map of each barcode ID and"
    echo " sequence is located in bdplex/index.192.txt"
    echo ""
    echo "Usage:"    
    echo "This script expects the FASTQ with all cells as the first argument, "
    echo " the FASTQ extension as the second e.g. _001.fastq.gz, .fq.gz;"
    echo " the output directory name as the third argument"
    echo ""
    echo "resource check w/ `/usr/bin/time -v`"
    echo "wannaAln:   2000 Kb == 2Mb"
    echo "zcat | awk | wc -l : 600 Kb .6 Mb"
    echo "bsub memory request set to 4 Gb"
    echo ""
    echo "Example:"
    echo "bsub -n 1 -M 4 -W 89 ./src/bdplex/barcode.wsplit.sh path/to/Sample_R1_001.fastq.gz _R1_001.fastq.gz wsplit/"
    echo "bsub -J 'wsplit[1-1500]' -n 1 -M 4 -W 356 ./src/bdplex/barcode.wsplit.sh path/to/Sample_R1_001.fastq.gz _R1_001.fastq.gz wsplit/"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 1 -M 4 -W 89 LSB_JOBINDEX=9 ./src/bdplex/barcode.wsplit.sh path/to/Sample_R1_001.fastq.gz _R1_001.fastq.gz wsplit/"
    echo ""
    exit 1;
}
#</usage>
set -e -x -o pipefail -u

## Argument 1 is file.R1.fastq.gz, 
## R2 is created by substituting _R1 with _R2
R1=$1
R2=$( echo ${R1} | sed -e 's/_R1/_R2/' )

## Argument 2 is the file extension to be replaced e.g. ".fastq.gz", ".fq.gz", "_IGO_XXXXX_A_001_R1.fq.gz"
EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fastq.gz

## Argument 3 is output directory
## if directory does not exist, its created
OUT=$3
[[ ! -z "$OUT" ]] || OUT=wsplit/
[ -d $OUT ] || mkdir -p $OUT

## file basename, removes extension
## used as prefix for all deplexing
X1=$( basename $R1 $EXTENSION )

## barcode names and sequences
## hardcoded and provided here for B-WGA Sets 1-96, and 97 to 192 in a single file
## only one barcode is searched at a time assigned byt the '$LSB_JOBINDEX' variable
BWGA192_IDX=( $( cut -f 1 src/bdeplex/index.192.txt ) )
BWGA192_SEQ=( $( cut -f 2 src/bdeplex/index.192.txt ) )

## output barcode id and sequence
echo "Searching for" ${BWGA192_IDX[$LSB_JOBINDEX]} ${BWGA192_SEQ[$LSB_JOBINDEX]}

## if nreads file exists: exit
[[ ! -f ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.nreads.txt ]] || exit 1

## PE demultiplexing -- core command
## use of wannaAln
## note `-q` is the correct flag for query pattern, not `-p` as in website
## -m 0 = perfect matching
## checkpoint if command runs succesfully, an .ok file is produced
wannaAln -a ${R1} -b ${R2} -x $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz \
	 -y $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz \
	 -q ${BWGA192_SEQ[$LSB_JOBINDEX]} -m 0
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} barcode split complete" > $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi

## Read counts per barcode, both R1 and R2 are counted
echo ${X1}.${BWGA192_IDX[$LSB_JOBINDEX]} \
     $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz | awk 'NR%4==1' | wc -l )  \
     $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz | awk 'NR%4==1' | wc -l ) | \
    tr ' ' "\t" > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.nreads.txt
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} read count complete" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi
