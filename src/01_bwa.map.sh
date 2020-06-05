#!/bin/bash
#BSUB -n 8 -R "rusage[mem=5]" -W 359 -R "span[hosts=1]"
## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -e -x -E -u -o pipefail

MAX_MEM_GB=38G

echo "LSB_JOBINDEX: " $LSB_JOBINDEX
echo "LSB_MAX_NUM_PROCESSORS: " $LSB_MAX_NUM_PROCESSORS

BWA_h37=/ifs/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
GFF_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Annotation/Genes/gencode.v19.annotation.gtf


R1=( $( find $1 -name "*R1*" | sort ) )
#R2=( $( find $1 -name "*R2*" | sort ) )
echo ${R1[@]}
echo ${R2[@]}

EXTENSION=$2 || true
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fastq.gz

OUT=$3 || true
[[ ! -z "$OUT" ]] || OUT=bwa_map/
[[ -d $OUT ]] || mkdir $OUT

MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//' ).bulk

echo $MID

[[ -d log/${MID} ]] || mkdir -p log/${MID}
[[ -d tmp ]] || mkdir -p tmp/


## quality control with fastp
## [[ -d fastp/ ]] || mkdir fastp
## fastp -i ${R1[$LSB_JOBINDEX]} -I ${R2[$LSB_JOBINDEX]} -o fastp/${MID}.R1.fq.gz -O fastp/${MID}.R2.fq.gz -p -h fastp/${MID}.html -j fastp/${MID}.json -R ${MID}  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").fastp.log
## if [ $? -eq 0 ] ; then echo "fastp succesfull" ; touch fastp/${MID}.fastp.ok ; else exit 1 ; fi

## alignment with BWA
## bwa mem -aM -t $LSB_MAX_NUM_PROCESSORS $BWA_h37 fastp/${MID}.R1.fq.gz fastp/${MID}.R2.fq.gz | samtools view -bS - > $OUT/${MID}.bam  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").bwa.log
## if [ $? -eq 0 ] ; then echo "alignment succesful" ; else exit 1; fi

## sort bam file
## sambamba sort -t $LSB_MAX_NUM_PROCESSORS -m $MAX_MEM_GB --tmpdir=tmp/ $OUT/${MID}.bam  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").sort.log
## if [ $? -eq 0 ] ; then touch ${OUT}/${MID}.bam.ok ; else exit 2 ; fi

## remove duplicate reads
## sambamba markdup --remove-duplicates -t $LSB_MAX_NUM_PROCESSORS --tmpdir=tmp/ $OUT/${MID}.sorted.bam $OUT/${MID}.dd.bam   2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").markdup.log
## if [ $? -eq 0 ] ; then rm $OUT/${MID}.sorted.bam{,.bai} ; touch ${OUT}/${MID}.dd.bam.ok ; else exit 3 ; fi


## additional metrics and QC
#for i in fastp metrics flagstats idxstats stats bamqc fastqc fastq_screen ; do [ -d $i ] || mkdir $i ; done

## Aligned reads quality assesment
#picard CollectMultipleMetrics I=$OUT/${MID}.dd.bam O=metrics/$MID R=$BWA_h37  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").picard.log
#if [ $? -eq 0 ] ; then echo "metrics succesfull" ; else exit 4 ; fi

#samtools flagstat -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.dd.bam > flagstats/${MID}.flagstat  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").flagstat.log
#if [ $? -eq 0 ] ; then echo "samtools flagstat succesfull" ; else exit 5 ; fi

#samtools idxstats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.dd.bam > idxstats/${MID}.idxstats  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").idxstats.log
#if [ $? -eq 0 ] ; then echo "samtools idxstats succesfull" ; else exit 6 ; fi

#samtools stats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.dd.bam > stats/${MID}.samtools.stats  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").samtools.log
#if [ $? -eq 0 ] ; then echo "samtools stats succesfull" ; else exit 7 ; fi

## source ~/opt/miniconda3/bin/activate  fastq-manip

#unset DISPLAY
#qualimap bamqc -nt $LSB_MAX_NUM_PROCESSORS -bam $OUT/${MID}.dd.bam -gd HUMAN -outdir bamqc/${MID}/ -outfile ${MID}.pdf -outformat "PDF:HTML" --java-mem-size=${MAX_MEM_GB}  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").bamqc.log
#if [ $? -eq 0 ] ; then echo "qualimap succesfull" ; else exit 8 ; fi

## FASTQ quality assesment
#fastqc -t $LSB_MAX_NUM_PROCESSORS -o fastqc/ ${R1[$LSB_JOBINDEX]} ${R2[$LSB_JOBINDEX]}  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").fastqc.log
#if [ $? -eq 0 ] ; then echo "fastqc succesfull" ; else exit 9 ; fi

fastq_screen --conf ~/conf/fastq_screen.conf --threads $LSB_MAX_NUM_PROCESSORS --top 8000000,100000 --outdir fastq_screen/ ${R1[$LSB_JOBINDEX]}  ${R2[$LSB_JOBINDEX]} --force  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").fastqscreen.log
if [ $? -eq 0 ] ; then echo "fastq_screen succesfull" ; else exit 7 ; fi

## fastq_screen --conf ~/conf/fastq_screen.conf --threads $LSB_MAX_NUM_PROCESSORS --top 3000000,100000 --outdir fastq_screen/ ${R1[$LSB_JOBINDEX]}  ${R2[$LSB_JOBINDEX]} --force  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").fastqscreen.log
fastp -i ${R1[$LSB_JOBINDEX]} -I ${R2[$LSB_JOBINDEX]} -w $LSB_MAX_NUM_PROCESSORS -p -h fastp/${MID}.html -j fastp/${MID}.json -R ${MID}  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").fastp.log
if [ $? -eq 0 ] ; then echo "fastp succesfull" ; else exit 10 ; fi

#__end__
