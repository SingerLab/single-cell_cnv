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
    echo "Example:"
    echo ""
    echo " bsub -Is -n 3 -M 3 -W 389 ./src/03_vbData.sh WD5816 bowtie SE"
    echo " bsub -Is -n 3 -M 3 -W 389 ./src/03_vbData.sh WD5816 bwa PE"
    echo ""
    exit 1;
}
#</usage>

## don't set -x file counts wont render properly
set -e -o pipefail -u

MID=${1}
ALIGNER=${2}
ALIGN_TYPE=${3}

[[ ! -z $MID ]] ||  { echo "Sample ID missing!" ; exit 1; }

echo "aggregating files:"
if [ ${ALIGN_TYPE} = "SE" ] ; then
    echo "# number of files in :"
    echo "#bin.dir vbData ploidy vbStats" | tr ' ' "\t"
    for i in  varbin{5,20,50}k ; do
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
	Rscript ./src/03_vbData.R  --sample.name=${MID} --input.dir=varbin${i} --output.dir=vbData --bin.size=${i} --aligner=${ALIGNER} &
    done

    echo "##" `date` >> processed.files.txt
    ./src/00_processed.files.sh $MID wsplit/ bowtie_out/ >> processed.files.txt
    wait ${!}
fi

if [ ${ALIGN_TYPE} = "PE" ] ; then
    echo "# number of files in :"
    echo "#bin.dir vbData ploidy vbStats" | tr ' ' "\t"
    for i in  varbin{5,20}k ; do
	echo ${i}/ $(ls $i/*bwa_seg_nobad.txt |wc -l ) \
	     $( ls $i/*quantal.ploidy.txt | wc -l ) \
	     $( ls $i/*.counts.stats.bed | wc -l ) ;
    done | tr ' ' "\t"

    BINS=(
	20k
	5k
    )
    for i in ${BINS[@]}
    do
	Rscript ./src/03_vbData_pe.R  --sample.name=${MID} --input.dir=varbin${i} --output.dir=vbData --bin.size=${i} --aligner=${ALIGNER} &
    done

    echo "##" `date` >> processed.files.txt
    ./src/00_processed.files.sh $MID wsplit_pe/ bwa_out/ >> processed.files.txt
    wait ${!}
fi

echo "sequence.quality per.sequence.quality per.sequence.gc per.base.N sequence.duplication overrepresented.sequences" | tr " " "\t" > vbData/${MID}_exclude.failed.txt
echo "sequence.quality per.sequence.quality per.sequence.gc per.base.N sequence.duplication overrepresented.sequences" | tr " " "\t" > vbData/${MID}_exclude.warnings.txt



## echo "starting geneCN"
## Rscript ./src/04_geneCN.R --sample.name=${MID} --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &
echo "estimating FGA"
Rscript ./src/04_fga.R --sample.name=${MID} --bin.size=5k --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &
## echo "simple heatmap"
## Rscript ./src/06_vbHeatmap.R --sample.name=${MID} --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &

wait ${!}

