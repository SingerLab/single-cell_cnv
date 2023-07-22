#!/bin/bash
#BSUB -n 3 -M 4 -W 359
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This is a wrapper script to prepare aggregated matrices containing"
    echo " integer copy number, lowess-ratio, and QC metrics"
    echo ""
    echo "Usage:"    
    echo "This script expects a common bioID as the first argument, and the aligner"
    echo " as a second argument.  The script expects the varbin{5,20,50}k/ directories"
    echo " from our single-cell_cnv pipeline"
    echo ""
    echo " bsub -Is -n 3 -M 3 -W 389 ./src/03_vbData.sh <genome> <bioID> <aligner> <SE|PE> <nobad>"
    echo ""
    echo " Where <genome> reference genome used, e.g. grch37, grch38, mmu10,"
    echo " <bioID> is a common biological ID across cells from the same sample/patient/project, etc,"
    echo " <aligne> aligner used"
    echo " <SE|PE> was sequencing single-end or paired-end"
    echo " <nobad> if TRUE, use the .*seg_nobad.txt file. else use the .*seg.txt"
    echo ""
    echo "Example:"
    echo ""
    echo " bsub -Is -n 3 -M 3 -W 389 ./src/03_vbData.sh grch38 WD5816 bowtie SE FALSE"
    echo " bsub -Is -n 3 -M 3 -W 389 ./src/03_vbData.sh grch38 WD5816 bwa PE FALSE"
    echo ""
    exit 1;
}
#</usage>

## don't set -x file counts wont render properly
set -e -o pipefail -u

GENOME=${1}
MID=${2}
ALIGNER=${3}
ALIGN_TYPE=${4}
NOBAD=${5}

[[ ! -z $MID ]] ||  { echo "Sample ID missing!" ; exit 1; }


if [ ${NOBAD} = "TRUE" ] ; then
    EXT=bwa_seg_nobad.txt
else 
    EXT=bwa_seg.txt
fi

echo "aggregating files:"
if [ ${ALIGN_TYPE} = "SE" ] ; then
    echo "# number of files in :"
    echo "#bin.dir vbData ploidy vbStats" | tr ' ' "\t"
    for i in  varbin{5,20,50,100}k ; do
	echo ${i}/ $(ls $i/*k50.varbin.data.txt |wc -l ) \
	     $( ls $i/*k50.varbin.quantal.stats.txt | wc -l ) \
	     $( ls $i/*.varbin.stats.txt | wc -l ) ;
    done | tr ' ' "\t"
    
    BINS=(
	50k
	20k
	5k )
    
    for i in ${BINS[@]}
    do
	Rscript ./src/03_vbData.R  --sample.name=${MID} --input.dir=varbin${i}/${MID} --output.dir=vbData/${MID} --genome=${GENOME} --bin.size=${i} --aligner=${ALIGNER} --nobad=${NOBAD} &
    done

    wait ${!}

    echo "##" `date` >> processed.files.txt
    ./src/00_processed.files.sh $MID wsplit/ bowtie_out/ >> processed.files.txt
fi

if [ ${ALIGN_TYPE} = "PE" ] ; then
    echo "# number of files in :"
    echo "#bin.dir vbData ploidy vbStats" | tr ' ' "\t"

    BINS=(
	100k
	50k
	20k
	5k
    )

    for i in ${BINS[@]} ; do
	echo varbin${i}/ $(ls varbin${i}/*${EXT} |wc -l ) \
	     $( ls varbin${i}/*quantal.ploidy.txt | wc -l ) \
	     $( ls varbin${i}/*.counts.stats.bed | wc -l ) ;
    done | tr ' ' "\t"
    
    for i in ${BINS[@]}
    do
	Rscript ./src/03_vbData_pe.R  --sample.name=${MID} --input.dir=varbin${i}/${MID} --genome=${GENOME} --nobad=${NOBAD} --output.dir=vbData/${MID} --bin.size=${i} --aligner=${ALIGNER} &
    done

    wait ${!}
    
    ## echo "##" `date` >> processed.files.txt
    ## ./src/00_processed.files.sh $MID wsplit_pe/ bwa_out/ >> processed.files.txt

fi

echo "sequence.quality per.sequence.quality per.sequence.gc per.base.N sequence.duplication overrepresented.sequences" | tr " " "\t" > vbData/${MID}/${MID}_exclude.failed.txt
echo "sequence.quality per.sequence.quality per.sequence.gc per.base.N sequence.duplication overrepresented.sequences" | tr " " "\t" > vbData/${MID}/${MID}_exclude.warnings.txt


## echo "starting geneCN"
## Rscript ./src/04_geneCN.R --sample.name=${MID} --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &
## echo "estimating FGA"
## Rscript ./src/04_fga.R --sample.name=${MID} --bin.size=5k --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &
## echo "simple heatmap"
## Rscript ./src/06_vbHeatmap.R --sample.name=${MID} --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &

## wait ${!}

