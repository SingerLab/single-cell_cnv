#!/bin/bash
#BSUB -n 3 -R "rusage[mem=3]" -R "hosts[span=1]" -W 358
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This script to runs the VARBIN algorithm on 5k, 20k, and 50k bins in the"
    echo " GRCh37 reference genome.  It calculates ploidy within the specified range"
    echo " using LOWESS least square means (ginkgo single-cell method)."
    echo ""
    echo "Usage:"    
    echo "This script expects a bam file as the first argument, file extension as the"
    echo " second argument; and the ploidy range as \`min' \`max' on the third and"
    echo " fourth arguments, respectively."
    echo ""
    echo "Example:"
    echo "bsub -n 3 -M 3 -W 89 ./src/02_varbin_pe.sh path/to/file.dd.bam .bam 1.5 4.8"
    echo ""
    exit 1;
}
#</usage>

## suggestions from:
## https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -e -x -o pipefail -u

cna_utils=$HOME/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/cna_utils/

BAM=$1

OUT5=varbin5k/
[ -d $OUT5 ] || mkdir -p $OUT5
OUT20=varbin20k/
[ -d $OUT20 ] || mkdir -p $OUT20
## OUT50=varbin50k
## [ -d $OUT50 ] || mkdir -p $OUT50

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.bam

MID=$( basename $BAM $EXTENSION )
echo $MID

## ploidy multipliers, min and max
MIN=$3
MAX=$4


$cna_utils/scripts/getBinCounts.py -i $BAM -b $cna_utils/data/hg19_5k_gz_enc_bins.bed \
    -d $cna_utils/data/hg19_150bp_dz_enc.bed -o ${OUT5}/${MID}_grch37.5k.bwa.bin.counts.bed \
    -v > ${OUT5}/${MID}_grch37.5k.bwa.bin.counts.stats.bed

$cna_utils/scripts/getBinCounts.py -i $BAM -b $cna_utils/data/hg19_20k_gz_enc_bins.bed \
    -d $cna_utils/data/hg19_150bp_dz_enc.bed -o ${OUT20}/${MID}_grch37.20k.bwa.bin.counts.bed \
    -v > ${OUT20}/${MID}_grch37.20k.bwa.bin.counts.stats.bed

wait


module load R/R-4.0.5

$cna_utils/scripts/cnvProfile.R -b ${OUT5}/${MID}_grch37.5k.bwa.bin.counts.bed -g $cna_utils/data/hg19_5k_gz_enc_gc.bed \
				-e $cna_utils/data/hg19_5k_gz_enc_badbins.bed \
				-n ${OUT5}/${MID}_grch37.5k.bwa -v \
				--minploidy=${MIN} --maxploidy=${MAX} 2> ${OUT5}/${MID}_grch37.5k.quantal.log

echo "cellID ploidy error" | tr ' ' '\t' >  ${OUT5}/${MID}_grch37.5k.quantal.ploidy.txt
echo $MID \
     $(grep "Ploidy" ${OUT5}/${MID}_grch37.5k.quantal.log  | tail -n 1 | sed -e 's/.*Ploidy: //') \
     $(grep "Error" ${OUT5}/${MID}_grch37.5k.quantal.log  | tail -n 1 | sed -e 's/.*Error: //') | \
    tr " " "\t" >> ${OUT5}/${MID}_grch37.5k.quantal.ploidy.txt

				
$cna_utils/scripts/cnvProfile.R -b ${OUT20}/${MID}_grch37.20k.bwa.bin.counts.bed -g $cna_utils/data/hg19_20k_gz_enc_gc.bed \
				-e $cna_utils/data/hg19_20k_gz_enc_badbins.bed \
				-n ${OUT20}/${MID}_grch37.20k.bwa -v \
				--minploidy=${MIN} --maxploidy=${MAX} 2> ${OUT20}/${MID}_grch37.20k.quantal.log

echo "cellID ploidy error" | tr ' ' '\t' >  ${OUT20}/${MID}_grch37.20k.quantal.ploidy.txt
echo $MID \
     $(grep "Ploidy" ${OUT20}/${MID}_grch37.20k.quantal.log  | tail -n 1 | sed -e 's/.*Ploidy: //') \
     $(grep "Error" ${OUT20}/${MID}_grch37.20k.quantal.log  | tail -n 1 | sed -e 's/.*Error: //') | \
    tr " " "\t" >> ${OUT20}/${MID}_grch37.20k.quantal.ploidy.txt
				

## __EOF__
