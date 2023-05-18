#!/bin/bash
#BSUB -n 4 -M 4 -W 89

#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This script to runs a second de-multiplexing step to generate pairwise "
    echo " paired-end (R1 and R2) reads for all possible barcode combinations "
    echo " of  each cell."
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
    echo "timing: 40 s by 1M reads"
    echo " estimate 1m per 1M reads for demultiplexing"
    echo ""
    echo "Example:"
    echo "bsub -n 1 -M 4 -W 89 ./src/bdplex/barcode.wsplit2X.sh path/to/cell.b99R1.fq.gz .R1.fq.gz wsplit2x/"
    echo "bsub -J 'wsplit[1-1500]' -n 1 -M 4 -W 356 ./src/bdplex/barcode.wsplit2X.sh path/to/cell.b99.R1.fq.gz .R1.fq.gz wsplit2x/"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 1 -M 4 -W 89 LSB_JOBINDEX=9 ./src/bdplex/barcode.wsplit2X.sh path/to/cell.b99.R1.fq.gz .R1.fq.gz wsplit2x/"
    echo ""
    exit 1;
}
#</usage>
set -e -x -o pipefail -u

## qc tool
FASTQC=/work/singer/opt/miniconda3/bin/fastqc
## genome reference
GRCH37="$HOME/cmo_genomes/GRCh37/bwa_fasta/b37.fasta"

## argument 1 R1 file
## read 1 input
## read 2 file created by substituting R1 to R2
R1=$1
R2=$( echo ${R1} | sed -e 's/.R1./.R2./' )

## extension specification argument 2
EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fq.gz

## X1 -- basename of cell 
X1=$( basename $R1 $EXTENSION )

## barcode 1 substring
BC1=$( echo $X1 | sed -e 's/.*b/b/' )

## barcode indices
BWGA192_IDX=( $( cut -f 1 src/bdeplex/index.192.txt ) )
BWGA192_SEQ=( $( cut -f 2 src/bdeplex/index.192.txt ) )

## barcode 2 substring
BC2=${BWGA192_IDX[$LSB_JOBINDEX]}

## echoing information
echo "Searching for" ${BWGA192_IDX[$LSB_JOBINDEX]} ${BWGA192_SEQ[$LSB_JOBINDEX]}

## output directory
OUT=$3
[[ ! -z "$OUT" ]] || OUT=wsplit2x/
[ -d $OUT ] || mkdir -p $OUT

## -- end of required argument specification

## checkpoint
## if nreads file exists: exit
[[ ! -f ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.chr.counts.txt ]] || exit 1

## build blank chromosome list to use when cell counts file is 0
[[ ! -f $(dirname $OUT)/grch37.chr.txt ]] || cut -f 1 ${GRCH37}.fai | sort -V > $(dirname $OUT)/grch37.chr.txt

## demultiplexing previously single-end demultiplexd cell against 1 barcode
##  specified by $LSB_JOBINDEX
## use of wannaAln
wannaAln -a ${R1} -b ${R2} -x $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz \
	 -y $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz \
	 -q ${BWGA192_SEQ[$LSB_JOBINDEX]} -m 0
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} barcode split complete" > $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi

## total file read count
echo ${X1}.${BWGA192_IDX[$LSB_JOBINDEX]} \
     $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz | awk 'NR%4==1' | wc -l )  \
     $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz | awk 'NR%4==1' | wc -l ) | \
    tr ' ' "\t" > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.nreads.txt
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} read count complete" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi

## number of reads
NR=$(cut -f 2 ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.nreads.txt )

## read counts per chromosome
## if read couts is zero, we set 0 to the list of all chromosomes
## only files with reads get processed, and each output file has only counts for
## the chromosomes prsent in the run.
## chunck streams the alignment and only counts chromosome reads
if [ $NR -gt 0 ] ; then     
    bwa mem -t ${LSB_MAX_NUM_PROCESSORS} -aM $GRCH37 $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz  | samtools view - | cut -f 3 | sort -V | uniq -c | awk '{$1=$1;print}' | tr ' ' '\t' > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.chr.counts.txt
else
    paste <(for i in {1..86} ; do echo 0 ; done) $(dirname ${OUT})/grch37.chr.txt > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.chr.counts.txt
fi
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} per chromosome counts complete" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi

## deleting unexpected barcode
## only perfect matching barcode files (expected pairing) get retained
if [ $BC1 != $BC2 ] ## && [ -f $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.chr.counts.txt ]
then
    rm $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz
    if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} R1 R2 fastq files deleted" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi
else
    echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} matching barcodes: R1 R2 fastq files kept" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok 
fi



