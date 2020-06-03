#!/bin/bash
#BSUB -n 8 -R "rusage[mem=4]" -R "span[hosts=1]" -W 359
## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -e -x -E -u -o pipefail

MAX_MEM_GB=32G

echo $LSB_JOBINDEX

echo $LSB_MAX_NUM_PROCESSORS

BWA_h37=/ifs/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
GFF_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Annotation/Genes/gencode.v19.annotation.gtf
BOWTIE_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/BowtieIndex/b37.fasta

R1=( $( find ${1} -name '*.gz' | grep R1 | sort) )

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fq.gz

OUT=$3
[[ ! -z "$OUT" ]] || OUT=bowtie_map/
[ -d $OUT ] || mkdir $OUT

MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//').bulk
echo $MID

[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/

## alignment with bowtie
bowtie --threads $LSB_MAX_NUM_PROCESSORS --sam --time --chunkmbs 512 -m 1 --best --strata $BOWTIE_h37 <( zcat ${R1[$LSB_JOBINDEX]} ) | samtools view -h -bS - > ${OUT}/${MID}.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").bowtie.log
if [ $? -eq 0 ] ; then echo "alignment succesful" ; touch ${OUT}/${MID}.bam.ok  ; else exit 1 ; fi

sambamba sort -m $MAX_MEM_GB --tmpdir=tmp/ -t $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").sambamba_sort.log
if [ $? -eq 0 ] ; then rm $OUT/${MID}.bam ; touch ${OUT}/${MID}.bam.ok ; else exit 2 ; fi

sambamba markdup --remove-duplicates --tmpdir=tmp/ -t $LSB_MAX_NUM_PROCESSORS --tmpdir=./ $OUT/${MID}.sorted.bam $OUT/${MID}.dd.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").sambamba_markdup.log
if [ $? -eq 0 ] ; then rm $OUT/${MID}.sorted.bam{,.bai} ${OUT}/${MID}.bam.ok ; touch ${OUT}/${MID}.dd.bam.ok ; else exit 3 ; fi

## submit to cn analysis with varbin
bsub -J "vb_$OUT/${MID}.dd.bam" src/02_varbin.sh $OUT/${MID}.dd.bam && echo "submitted to varbin" > log/${MID}/$(date "+%Y%m%d-%H%M%S").varbin_submit.log 2>&1

## additional metrics and QC
for i in fastp metrics idxstats stats fastq_screen flagstats preseq bamqc fastqc ; do [ -d $i ] || mkdir $i ; done

samtools flagstat -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.dd.bam > flagstats/${MID}.bowtie.flagstat  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_flagstat.log
if [ $? -eq 0 ] ; then echo "samtools flagstat succesfull" ; else exit 5 ; fi

samtools idxstats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.dd.bam > idxstats/${MID}.bowtie.idxstats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log
if [ $? -eq 0 ] ; then echo "samtools idxstats succesfull" ; else exit 6 ; fi

samtools stats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.dd.bam > stats/${MID}.bowtie.samtools.stats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_stats.log
if [ $? -eq 0 ] ; then echo "samtools stats succesfull" ; else exit 7 ; fi

unset DISPLAY
qualimap bamqc -nt $LSB_MAX_NUM_PROCESSORS -bam $OUT/${MID}.dd.bam -gd HUMAN -outdir bamqc/${MID}/ -outfile ${MID}.bowtie.pdf -outformat "PDF:HTML"  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").qualimap_bamqc.log
if [ $? -eq 0 ] ; then echo "qualimap succesfull" ; else exit 8 ; fi

preseq c_curve -B $OUT/${MID}.dd.bam > preseq/${MID}.preseq.c_curve 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.c_curve.log
if [ $? -eq 0 ] ; then echo "preseq c_curve succesfull" ; else exit 9 ; fi

preseq lc_extrap -B $OUT/${MID}.dd.bam > preseq/${MID}.preseq.lc_extrap 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.lc_extrap.log
if [ $? -eq 0 ] ; then echo "preseq c_curve succesfull" ; else exit 10 ; fi

## fastq quality assesment
fastp -i ${R1[$LSB_JOBINDEX]} -p -h fastp/${MID}.html -j fastp/${MID}.json -R ${MID} --umi --umi_loc read1 --umi_len 6  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastp.log
if [ $? -eq 0 ] ; then echo "fastp succesfull" ; touch fastp/${MID}.fastp.ok ; else exit 11 ; fi

fastqc -t $LSB_MAX_NUM_PROCESSORS -a ~/genomes/adapters/bwga/adapter_list.txt -o fastqc/ ${R1[$LSB_JOBINDEX]}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastqc.log
if [ $? -eq 0 ] ; then echo "fastqc succesfull" ; else exit 12 ; fi

fastq_screen --aligner bwa --threads $LSB_MAX_NUM_PROCESSORS --subset 0 --conf ~/conf/fastq_screen.conf --outdir fastq_screen/ ${R1[$LSB_JOBINDEX]} --force  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastq_screen.log
if [ $? -eq 0 ] ; then echo "fastqc_screen succesfull" ; else exit 13 ; fi

picard CollectMultipleMetrics I=$OUT/${MID}.dd.bam O=metrics/$MID R=${BOWTIE_h37}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").picard_metrics.log
if [ $? -eq 0 ] ; then echo "picard metrics succesfull" ; else exit 4 ; fi


#__end__#
