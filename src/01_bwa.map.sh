#!/bin/bash
#BSUB -n 8 -R "rusage[mem=5]" -W 359 -R "span[hosts=1]"
## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -e -x -E -u -o pipefail

MAX_MEM_GB=38G

echo "LSB_JOBINDEX: " $LSB_JOBINDEX
echo "LSB_MAX_NUM_PROCESSORS: " $LSB_MAX_NUM_PROCESSORS

BWA_h37=/ifs/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
GFF_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Annotation/Genes/gencode.v19.annotation.gtf


R1=( $( find -L $1 -name "*R1*" | sort ) )
R2=( $( find -L $1 -name "*R2*" | sort ) )
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

echo $MID

[[ -d log/${MID} ]] || mkdir -p log/${MID}
[[ -d tmp ]] || mkdir -p tmp/


## quality control with fastp
[[ -d fastp/ ]] || mkdir fastp
fastp -i ${R1[$LSB_JOBINDEX]} -I ${R2[$LSB_JOBINDEX]} -o fastp/${MID}.R1.fastq.gz -O fastp/${MID}.R2.fastq.gz -p -h fastp/${MID}.html -j fastp/${MID}.json -R ${MID}  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").fastp.log
if [ $? -eq 0 ] ; then echo "fastp succesfull" ; touch fastp/${MID}.fastp.ok ; else exit 1 ; fi

## alignment with BWA
bwa mem -aM -t $LSB_MAX_NUM_PROCESSORS $BWA_h37 fastp/${MID}.R1.fastq.gz fastp/${MID}.R2.fastq.gz | samtools view -bS - > $OUT/${MID}.bam  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").bwa.log
if [ $? -eq 0 ] ; then echo "alignment succesful" ; else exit 1; fi

## sort bam file
sambamba sort -t $LSB_MAX_NUM_PROCESSORS -m $MAX_MEM_GB --tmpdir=tmp/ $OUT/${MID}.bam  2> log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").sort.log
if [ $? -eq 0 ] ; then touch ${OUT}/${MID}.bam.ok ; else exit 2 ; fi

[[ -d metrics ]] || mkdir -p metrics/
## remove duplicate reads
picard MarkDuplicates I=$OUT/${MID}.sorted.bam O=$OUT/${MID}.md.bam M=metrics/${MID}.picard.markdups  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").picard_markdup.log
if [ $? -eq 0 ] ; then rm $OUT/${MID}.sorted.bam{,.bai} ${OUT}/${MID}.bam{,.ok} ; else exit 6 ; fi

samtools index -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.md.bam




#__end__
