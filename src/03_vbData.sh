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
    echo " bsub -Is -n 3 -M 3 -W 389 ./src/03_vbData.sh WD5816 bowtie"
    echo ""
    exit 1;
}
#</usage>

echo "aggregating files:"
echo "# number of files in :"
echo "#bin.dir vbData ploidy vbStats" | tr ' ' "\t"
for i in  varbin{5,20,50}k ; do echo ${i}/ $(ls $i/*k50.varbin.data.txt |wc -l ) $( ls $i/*k50.varbin.quantal.stats.txt | wc -l ) $( ls $i/*.varbin.stats.txt | wc -l ) ; done | tr ' ' "\t"

set -e -x -o pipefail -u

MID=${1}
ALIGNER=${2}

[[ ! -z $MID ]] ||  { echo "Sample ID missing!" ; exit 1; }

for i in 50k 20k 5k
do
    Rscript ./src/03_vbData.R  --sample.name=${MID} --input.dir=varbin${i} --output.dir=vbData --bin.size=${i} --aligner=${ALIGNER} &
done

echo "sequence.quality per.sequence.quality per.sequence.gc per.base.N sequence.duplication overrepresented.sequences" | tr " " "\t" > vbData/${MID}_exclude.failed.txt
echo "sequence.quality per.sequence.quality per.sequence.gc per.base.N sequence.duplication overrepresented.sequences" | tr " " "\t" > vbData/${MID}_exclude.warnings.txt

wait ${!}

## echo "starting geneCN"
## Rscript ./src/04_geneCN.R --sample.name=${MID} --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &
echo "estimating FGA"
Rscript ./src/04_fga.R --sample.name=${MID} --bin.size=5k --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &
echo "simple heatmap"
Rscript ./src/06_vbHeatmap.R --sample.name=${MID} --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &

wait ${!}

set +x

echo "##" `date` >> processed.files.txt
./src/00_processed.files.sh $MID wsplit/ bowtie_out/ >> processed.files.txt
