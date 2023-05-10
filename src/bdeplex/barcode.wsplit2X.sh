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
    echo "bsub -n 1 -M 4 -W 89 ./src/bdplex/barcode.wsplit2X.sh path/to/cell.b99R1.fq.gz .R1.fq.gz wsplit_2x/"
    echo "bsub -J 'wsplit[1-1500]' -n 1 -M 4 -W 356 ./src/bdplex/barcode.wsplit2X.sh path/to/cell.b99.R1.fq.gz .R1.fq.gz wsplit_2x/"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 1 -M 4 -W 89 LSB_JOBINDEX=9 ./src/bdplex/barcode.wsplit2X.sh path/to/cell.b99.R1.fq.gz .R1.fq.gz wsplit_2x/"
    echo ""
    exit 1;
}
#</usage>
set -e -x -o pipefail -u

FASTQC=/work/singer/opt/miniconda3/bin/fastqc
GRCH37="$HOME/cmo_genomes/GRCh37/bwa_fasta/b37.fasta"

R1=$1
R2=$( echo ${R1} | sed -e 's/.R1./.R2./' )

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fq.gz

X1=$( basename $R1 $EXTENSION )
BC1=$( echo $X1 | sed -e 's/.*b/b/' )

BWGA192_IDX=( $( cut -f 1 src/bdeplex/index.192.txt ) )
BWGA192_SEQ=( $( cut -f 2 src/bdeplex/index.192.txt ) )

BC2=${BWGA192_IDX[$LSB_JOBINDEX]}

echo "Searching for" ${BWGA192_IDX[$LSB_JOBINDEX]} ${BWGA192_SEQ[$LSB_JOBINDEX]}

OUT=$3
[[ ! -z "$OUT" ]] || OUT=wsplit_2x/
[ -d $OUT ] || mkdir -p $OUT


## if nreads file exists: exit
[[ ! -f ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.chr.counts.txt ]] || exit 1

## build blank reference in case counts file is 0
[[ ! -f $(dirname $OUT)/grch37.chr.txt ]] || cut -f 1 ${GRCH37}.fai | sort -V > $(dirname $OUT)/grch37.chr.txt

## use of wannaAln
wannaAln -a ${R1} -b ${R2} -x $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz \
	 -y $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz \
	 -q ${BWGA192_SEQ[$LSB_JOBINDEX]} -m 0
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} barcode split complete" > $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi

echo ${X1}.${BWGA192_IDX[$LSB_JOBINDEX]} \
     $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz | awk 'NR%4==1' | wc -l )  \
     $( zcat ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz | awk 'NR%4==1' | wc -l ) | \
    tr ' ' "\t" >> ${OUT}/cell.nreads.txt
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} read count complete" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi

NR=$(cut -f 2 ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.nreads.txt )

if [ $NR -gt 0 ] ; then     
    bwa mem -t ${LSB_MAX_NUM_PROCESSORS} -aM $GRCH37 $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz  | samtools view - | cut -f 3 | sort -V | uniq -c | awk '{$1=$1;print}' | tr ' ' '\t' > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.chr.counts.txt
else
    paste <(for i in {1..86} ; do echo 0 ; done) $(dirname ${OUT})/grch37.chr.txt > ${OUT}/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.chr.counts.txt
fi
if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} per chromosome counts complete" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi


if [ $BC1 != $BC2 ] ## && [ -f $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.chr.counts.txt ]
then
    rm $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R1.fq.gz $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.R2.fq.gz
    if [ $? -eq 0 ] ; then echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} R1 R2 fastq files deleted" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok ; fi
else
    echo `date` ": ${BWGA192_IDX[$LSB_JOBINDEX]} matching barcodes: R1 R2 fastq files kept" >> $OUT/${X1}.${BWGA192_IDX[$LSB_JOBINDEX]}.ok 
fi





########################################################################
## -- old script
########################################################################
#% #<usage>
#% [[ $# -gt 0 ]] || {
#%     echo "Description:"
#%     echo "Identifying corcodant barcodes in R1 and R2 files"
#%     echo "This script requires to run the barcode.wsplit.sh first"
#%     exit 1;
#% }
#% #</usage>
#% set -e -x -o pipefail -u
#% 
#% R1=$1
#% R2=$( echo ${R1} | sed -e 's/_R1/_R2/' )
#% 
#% EXTENSION=$2
#% [[ ! -z "$EXTENSION" ]] || EXTENSION=.fastq.gz
#% 
#% X1=$( basename $R1 $EXTENSION )
#% 
#% BC=$( echo $X1 | sed -e 's/.*b//' )
#% echo "Barcode is" $BC
#% 
#% if [ $BC -eq 1 ]
#% then
#%     BCIDX=0
#% elif [ $BC -ne 1 ]
#% then
#%     BCIDX=$(expr $BC - 1)
#% fi
#% 
#% echo $BCIDX
#% 
#% BWGA192_IDX=( $( cut -f 1 src/bdeplex/index.192.txt ) )
#% BWGA192_SEQ=( $( cut -f 2 src/bdeplex/index.192.txt ) )
#% 
#% echo "Searching for" ${BWGA192_IDX[$BCIDX]} ${BWGA192_SEQ[$BCIDX]}
#% 
#% OUT=$3
#% [[ ! -z "$OUT" ]] || OUT=wsplit2/
#% [ -d $OUT ] || mkdir -p $OUT/
#% 
#% tmp_out_file=${OUT}/tmp.${X1}.${BWGA192_IDX[$BCIDX]}.R1.R2.txt.gz
#% 
#% R1_OUT_FILE=${OUT}/${X1}.${BWGA192_IDX[$BCIDX]}.R1.fq.gz
#% R2_OUT_FILE=${OUT}/${X1}.${BWGA192_IDX[$BCIDX]}.R2.fq.gz
#% 
#% paste <( zcat ${R1} | paste  - - - - )  \
#%       <( zcat ${R2} | paste  - - - - ) | \
#%     awk '$2 ~ /'${BWGA192_SEQ[$BCIDX]}'[ACTGN]{83}/ &&
#%          $6 ~ /'${BWGA192_SEQ[$BCIDX]}'[ACTGN]{83}/ {print}' | \
#%     gzip -c > ${tmp_out_file}
#% 
#% zcat ${tmp_out_file} | cut -f 1-4 | tr "\t" "\n" | gzip -c > $R1_OUT_FILE
#% zcat ${tmp_out_file} | cut -f 5-8 | tr "\t" "\n" | gzip -c > $R2_OUT_FILE
#% 
#% rm ${tmp_out_file}
#% 
#% echo ${X1}.${BWGA192_IDX[$BCIDX]} \
#%      $( zcat ${R1_OUT_FILE} | awk 'NR%4==1' | wc -l )  \
#%      $( zcat ${R2_OUT_FILE} | awk 'NR%4==1' | wc -l ) | \
#%     tr ' ' "\t" > ${OUT}/${X1}.${BWGA192_IDX[$BCIDX]}.nreads.txt
#% 
