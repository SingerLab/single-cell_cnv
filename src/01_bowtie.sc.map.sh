#!/bin/bash
#BSUB -n 8 -R "rusage[mem=4]" -R "span[hosts=1]" -W 89

## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -E -e -x -u -o pipefail

## trap read debug
MAX_MEM_GB=22G

echo $LSB_JOBINDEX
echo $LSB_MAX_NUM_PROCESSORS

BWA_h37=/ifs/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
GFF_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Annotation/Genes/gencode.v19.annotation.gtf
BOWTIE_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/BowtieIndex/b37.fasta

R1=( $( ls ${1}*fq.gz ) )

EXTENSION=$2
[[ ! -z "$EXTENSION" ]] || EXTENSION=.fq.gz

OUT=$3
[[ ! -z "$OUT" ]] || OUT=bowtie_map/
[ -d $OUT ] || mkdir $OUT

MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//')
echo $MID

MIN_PLOIDY=$4
[[ ! -z "$MIN_PLOIDY" ]] || MIN_PLOIDY=1.5

MAX_PLOIDY=$5
[[ ! -z "$MAX_PLOIDY" ]] || MAX_PLOIDY=4.5

[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/
[ -d metrics/ ] || mkdir -p metrics/

## if bam file exists: exit
[[ ! -f  $OUT/${MID}.md.bam ]] || exit 1

## 3 prime trim based on read length
## estimating read length
bpcount=$( echo $(zcat ${R1[$LSB_JOBINDEX]} | head -n 4 | awk 'NR==2 {print length}') )

## setting trim3
if [ $bpcount -le 51 ]
then
    trim5=0
    trim3=0
elif [ $bpcount -le 101 ]
then
    trim5=48
    trim3=0
elif [ $bpcount -eq 151 ]
then
    trim5=48
    trim3=50
fi

## alignment with bowtie
bowtie --threads $LSB_MAX_NUM_PROCESSORS --sam --time --chunkmbs 256 --trim5 $trim5 --trim3 $trim3 -m 1 --best --strata $BOWTIE_h37 <( zcat ${R1[$LSB_JOBINDEX]} ) | samtools view -h -bS - > ${OUT}/${MID}.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").bowtie.log
if [ $? -eq 0 ] ; then echo "alignment succesful" ; touch ${OUT}/${MID}.bam.ok  ; else exit 3 ; fi

sambamba sort -m $MAX_MEM_GB --tmpdir=tmp/ -t $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").sambamba_sort.log
if [ $? -eq 0 ] ; then rm $OUT/${MID}.bam ; touch ${OUT}/${MID}.bam.ok ; else exit 4 ; fi

## sambamba markdup --remove-duplicates --tmpdir=tmp/ -t $LSB_MAX_NUM_PROCESSORS --tmpdir=tmp/ $OUT/${MID}.sorted.bam $OUT/${MID}.dd.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").sambamba_markdup.log
picard MarkDuplicates I=$OUT/${MID}.sorted.bam O=$OUT/${MID}.dd.bam M=metrics/${MID}.picard.rmdups REMOVE_DUPLICATES=true  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").rmdup.log
if [ $? -eq 0 ] ; then touch ${OUT}/${MID}.dd.bam.ok ; else exit 5 ; fi
## indexing
samtools index $OUT/${MID}.dd.bam
if [ $? -eq 0 ] ; then touch ${OUT}/${MID}.dd.bam.ok ; else exit 5 ; fi

## submit to cn analysis with varbin
bsub -J "vb_$OUT/${MID}.dd.bam" src/02_varbin.sh $OUT/${MID}.dd.bam .bam $MIN_PLOIDY $MAX_PLOIDY > log/${MID}/$(date "+%Y%m%d-%H%M%S").varbin_submit.log 2>&1  && echo "submitted to varbin" 

## additional metrics and QC
for i in fastp metrics idxstats stats fastq_screen flagstats preseq bamqc fastqc ; do [ -d $i ] || mkdir $i ; done

picard MarkDuplicates I=$OUT/${MID}.sorted.bam O=$OUT/${MID}.md.bam M=metrics/${MID}.picard.markdups  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").picard_markdup.log
if [ $? -eq 0 ] ; then rm $OUT/${MID}.sorted.bam{,.bai} ${OUT}/${MID}.bam.ok ; else exit 6 ; fi

samtools flagstat -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.md.bam > flagstats/${MID}.bowtie.flagstat  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_flagstat.log
if [ $? -eq 0 ] ; then echo "samtools flagstat succesfull" ; else exit 7 ; fi

samtools idxstats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.md.bam > idxstats/${MID}.bowtie.idxstats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log
if [ $? -eq 0 ] ; then echo "samtools idxstats succesfull" ; else exit 8 ; fi

samtools stats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.md.bam > stats/${MID}.bowtie.samtools.stats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_stats.log
if [ $? -eq 0 ] ; then echo "samtools stats succesfull" ; else exit 9 ; fi

unset DISPLAY
qualimap bamqc -nt $LSB_MAX_NUM_PROCESSORS -bam $OUT/${MID}.md.bam -gd HUMAN -outdir bamqc/${MID}/ -outfile ${MID}.bowtie.pdf -outformat "PDF:HTML"  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").qualimap_bamqc.log
if [ $? -eq 0 ] ; then echo "qualimap succesfull" ; else exit 10 ; fi

preseq c_curve -B $OUT/${MID}.md.bam > preseq/${MID}.preseq.c_curve 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.c_curve.log
if [ $? -eq 0 ] ; then echo "preseq c_curve succesfull" ; else exit 11 ; fi

preseq lc_extrap -B $OUT/${MID}.md.bam > preseq/${MID}.preseq.lc_extrap 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.lc_extrap.log
if [ $? -eq 0 ] ; then echo "preseq c_curve succesfull" ; else exit 12 ; fi

## placed here as it often fails
picard CollectMultipleMetrics I=$OUT/${MID}.md.bam O=metrics/$MID R=${BOWTIE_h37}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").picard_metrics.log
if [ $? -eq 0 ] ; then echo "picard metrics succesfull" ; else exit 13 ; fi


#__end__#

## fastq quality assesment -- MOVED TO DIFFERENT FILES
## fastp -i ${R1[$LSB_JOBINDEX]} -p -h fastp/${MID}.html -j fastp/${MID}.json -R ${MID} --umi --umi_loc read1 --umi_len 6  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastp.log
## if [ $? -eq 0 ] ; then echo "fastp succesfull" ; touch fastp/${MID}.fastp.ok ; else exit 11 ; fi
## fastqc -t $LSB_MAX_NUM_PROCESSORS -a ~/genomes/adapters/bwga/adapter_list.txt -o fastqc/ ${R1[$LSB_JOBINDEX]}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastqc.log
## if [ $? -eq 0 ] ; then echo "fastqc succesfull" ; else exit 12 ; fi
## fastq_screen --aligner bwa --threads $LSB_MAX_NUM_PROCESSORS --top 1000000,1000 --conf ~/conf/fastq_screen.conf --outdir fastq_screen/ ${R1[$LSB_JOBINDEX]} --force  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastq_screen.log
## if [ $? -eq 0 ] ; then echo "fastqc_screen succesfull" ; else exit 13 ; fi

## for loops for failed QC 
## for i in bsplit/*gz ; do MID=$(basename $i .fq.gz ); bsub -n 8 -R "rusage[mem=6]" -W 89 fastq_screen --aligner bwa --threads 8 --subset 0 --conf ~/conf/fastq_screen.conf --outdir fastq_screen/ $i --force  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastq_screen.log ; done

## for i in fastq_screen/*fastq ; do ; in=bsplit/$( basename $i _temp_subset.fastq ) ; mid=$( basename $inFile .fq.gz ) ; bsub -n 8 -R "rusage[mem=6]" -W 89 fastq_screen --aligner bwa --threads 8 --subset 0 --conf ~/conf/fastq_screen.conf --outdir fastq_screen/ $i --force 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastq_screen.log ; done

## for i in bsplit/*gz ; do MID=$(basename $i .fq.gz ); bsub -n 8 -R "rusage[mem=6]" -W 89 fastq_screen --aligner bwa --threads 8 --subset 0 --conf ~/conf/fastq_screen.conf --outdir fastq_screen/ $i --force  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").fastq_screen.log ; done


## for i in bsplit/*gz ; do MID=$(basename $i .fq.gz ); bsub -n 1 -R "rusage[mem=12]" -W 24:00 "picard CollectMultipleMetrics I=bowtie_out/${MID}.md.bam O=metrics/$MID R=${BOWTIE_h37}  2> log/${MID}/$(date '+%Y%m%d-%H%M%S').picard_metrics.log" ; done

## for i in bowtie_out/*md.bam ; do mid=$( basename $i .md.bam) ; bsub -n 1 -M 8 -W 24:00 "picard CollectMultipleMetrics I=${i} O=metrics/$mid R=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/BowtieIndex/b37.fasta  2> log/${mid}/$(date '+%Y%m%d-%H%M%S').picard_metrics.log" ; done

## for i in $( cat metrics.missing.txt ); do MID=$i; bsub -n 1 -R "rusage[mem=12]" -W 24:00 picard CollectMultipleMetrics I=bowtie_out/${MID}.md.bam O=metrics/$MID R=${BOWTIE_h37} ; done


## for i in bsplit/*gz ; do mid=$( basename $i .fq.gz ) ; bsub -n 1 -M 4 -W 89 "fastp -i ${i} -p -h fastp/${mid}.html -j fastp/${mid}.json -R ${mid} --umi --umi_loc read1 --umi_len 6  2> log/${MID}/$(date '+%Y%m%d-%H%M%S').fastp.log" ; done

## for i in bsplit/*gz ; do mid=$( basename $i .fq.gz ) ; bsub -n 1 -M 4 -W 89 fastqc -t $LSB_MAX_NUM_PROCESSORS -a ~/genomes/adapters/bwga/adapter_list.txt -o fastqc/ ${i} 2> log/${mid}/$(date "+%Y%m%d-%H%M%S").fastqc.log ; done


