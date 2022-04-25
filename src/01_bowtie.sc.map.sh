#!/bin/bash
#BSUB -n 8 -R "rusage[mem=4]" -R "span[hosts=1]" -W 89
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This script to runs the base pipeline to align, FASTQ files to the human"
    echo " GRCh37 reference genome.  It then runs mark duplicates and duplicate"
    echo " read removal.  Furthermore, it runs QC metrics from picard-tools, and"
    echo " samtools on the bam files"
    echo ""
    echo "Usage:"    
    echo "This script expects a directory name with the FASTQ files as the first"
    echo " argument; the FASTQ extension as the second e.g. _001.fastq.gz, .fq.gz;"
    echo " the output directory name as the third argument; and the ploidy range"
    echo " as \`min' \`max' on the fourth and fifth arguments, respectively."
    echo ""
    echo "Example:"
    echo "bsub -n 8 -M 4 -W 89 ./src/01_bowtie.sc.map.sh path/to/fastq/ .fq.gz bowtie_out/ 1.5 4.8"
    echo "bsub -J 'bm[1-1500]' -n 8 -M 4 -W 89 ./src/01_bowtie.sc.map.sh path/to/fastq/ .fq.gz bowtie_out/ 1.5 4.8"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 8 -M 4 -W 89 LSB_JOBINDEX=9 ./src/01_bowtie.sc.map.sh path/to/fastq/ .fq.gz bowtie_out/ 1.5 4.8"
    echo ""
    exit 1;
}
#</usage>

## suggestions from https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -E -e -x -u -o pipefail

## trap read debug
MAX_MEM_GB=22G

echo $LSB_JOBINDEX
echo $LSB_MAX_NUM_PROCESSORS

BWA_h37=/juno/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
GFF_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Annotation/Genes/gencode.v19.annotation.gtf
BOWTIE_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/BowtieIndex/b37.fasta

## hsa.SE.dd.bam -- change to choose genomes
## DEDUP_BAM_EXT=hsa.SE.dd.bam
## MARKDUP_BAM_EXT=hsa.SE.md.bam

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
[[ ! -f  $OUT/${MID}.hsa.SE.md.bam ]] || exit 1

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
BOWTIE=/work/singer/opt/miniconda3/bin/bowtie
SAMTOOLS=/work/singer/opt/miniconda3/bin/samtools
SAMBAMBA=/home/gularter/bin/sambamba
PICARD=/work/singer/opt/miniconda3/bin/picard 
QUALIMAP=/work/singer/opt/miniconda3/bin/qualimap
PRESEQ=/home/gularter/bin/preseq

$BOWTIE --threads $LSB_MAX_NUM_PROCESSORS --sam --time --chunkmbs 256 --trim5 $trim5 --trim3 $trim3 -m 1 --best --strata $BOWTIE_h37 <( zcat ${R1[$LSB_JOBINDEX]} ) | $SAMTOOLS view -h -bS - > ${OUT}/${MID}.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").bowtie.log
if [ $? -eq 0 ] ; then echo "alignment succesful" ; touch ${OUT}/${MID}.bam.ok  ; else exit 3 ; fi

$SAMBAMBA sort -m $MAX_MEM_GB --tmpdir=tmp/ -t $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").sambamba_sort.log
if [ $? -eq 0 ] ; then rm $OUT/${MID}.bam ; touch ${OUT}/${MID}.bam.ok ; else exit 4 ; fi

$PICARD MarkDuplicates I=$OUT/${MID}.sorted.bam O=$OUT/${MID}.hsa.SE.dd.bam M=metrics/${MID}.picard.rmdups REMOVE_DUPLICATES=true  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").rmdup.log
if [ $? -eq 0 ] ; then touch ${OUT}/${MID}.hsa.SE.dd.bam.ok ; else exit 5 ; fi

## indexing
$SAMTOOLS index $OUT/${MID}.hsa.SE.dd.bam
if [ $? -eq 0 ] ; then touch ${OUT}/${MID}.hsa.SE.dd.bam.ok ; else exit 5 ; fi

## submit to cn analysis with varbin
bsub -J "vb_$OUT/${MID}.hsa.SE.dd.bam" src/02_varbin.sh $OUT/${MID}.hsa.SE.dd.bam .bam $MIN_PLOIDY $MAX_PLOIDY > log/${MID}/$(date "+%Y%m%d-%H%M%S").varbin_submit.log 2>&1  && echo "submitted to varbin" 

## additional metrics and QC
for i in fastp metrics idxstats stats fastq_screen flagstats preseq bamqc fastqc ; do [ -d $i ] || mkdir $i ; done

$PICARD MarkDuplicates I=$OUT/${MID}.sorted.bam O=$OUT/${MID}.hsa.SE.md.bam M=metrics/${MID}.picard.markdups  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").picard_markdup.log
if [ $? -eq 0 ] ; then rm $OUT/${MID}.sorted.bam{,.bai} ${OUT}/${MID}.bam.ok ; else exit 6 ; fi

$SAMTOOLS flagstat -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.hsa.SE.md.bam > flagstats/${MID}.bowtie.flagstat  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_flagstat.log
if [ $? -eq 0 ] ; then echo "samtools flagstat succesfull" ; else exit 7 ; fi

$SAMTOOLS idxstats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.hsa.SE.md.bam > idxstats/${MID}.bowtie.idxstats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log
if [ $? -eq 0 ] ; then echo "samtools idxstats succesfull" ; else exit 8 ; fi

$SAMTOOLS stats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.hsa.SE.md.bam > stats/${MID}.bowtie.samtools.stats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_stats.log
if [ $? -eq 0 ] ; then echo "samtools stats succesfull" ; else exit 9 ; fi

unset DISPLAY
$QUALIMAP bamqc -nt $LSB_MAX_NUM_PROCESSORS -bam $OUT/${MID}.hsa.SE.md.bam -gd HUMAN -outdir bamqc/${MID}/ -outfile ${MID}.bowtie.pdf -outformat "PDF:HTML"  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").qualimap_bamqc.log
if [ $? -eq 0 ] ; then echo "qualimap succesfull" ; else exit 10 ; fi

$PRESEQ c_curve -B $OUT/${MID}.hsa.SE.md.bam > preseq/${MID}.preseq.c_curve 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.c_curve.log
if [ $? -eq 0 ] ; then echo "preseq c_curve succesfull" ; else exit 11 ; fi

$PRESEQ lc_extrap -B $OUT/${MID}.hsa.SE.md.bam > preseq/${MID}.preseq.lc_extrap 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.lc_extrap.log
if [ $? -eq 0 ] ; then echo "preseq c_curve succesfull" ; else exit 12 ; fi

## placed here as it often fails
$PICARD CollectMultipleMetrics I=$OUT/${MID}.hsa.SE.md.bam O=metrics/$MID R=${BOWTIE_h37}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").picard_metrics.log
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


## for i in bsplit/*gz ; do MID=$(basename $i .fq.gz ); bsub -n 1 -R "rusage[mem=12]" -W 24:00 "picard CollectMultipleMetrics I=bowtie_out/${MID}.hsa.SE.md.bam O=metrics/$MID R=${BOWTIE_h37}  2> log/${MID}/$(date '+%Y%m%d-%H%M%S').picard_metrics.log" ; done

## for i in bowtie_out/*md.bam ; do mid=$( basename $i .hsa.SE.md.bam) ; bsub -n 1 -M 8 -W 24:00 "picard CollectMultipleMetrics I=${i} O=metrics/$mid R=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/BowtieIndex/b37.fasta  2> log/${mid}/$(date '+%Y%m%d-%H%M%S').picard_metrics.log" ; done

## for i in $( cat metrics.missing.txt ); do MID=$i; bsub -n 1 -R "rusage[mem=12]" -W 24:00 picard CollectMultipleMetrics I=bowtie_out/${MID}.hsa.SE.md.bam O=metrics/$MID R=${BOWTIE_h37} ; done


## for i in bsplit/*gz ; do mid=$( basename $i .fq.gz ) ; bsub -n 1 -M 4 -W 89 "fastp -i ${i} -p -h fastp/${mid}.html -j fastp/${mid}.json -R ${mid} --umi --umi_loc read1 --umi_len 6  2> log/${MID}/$(date '+%Y%m%d-%H%M%S').fastp.log" ; done

## for i in bsplit/*gz ; do mid=$( basename $i .fq.gz ) ; bsub -n 1 -M 4 -W 89 fastqc -t $LSB_MAX_NUM_PROCESSORS -a ~/genomes/adapters/bwga/adapter_list.txt -o fastqc/ ${i} 2> log/${mid}/$(date "+%Y%m%d-%H%M%S").fastqc.log ; done

## previous deduplication of reads
## sambamba markdup --remove-duplicates --tmpdir=tmp/ -t $LSB_MAX_NUM_PROCESSORS --tmpdir=tmp/ $OUT/${MID}.sorted.bam $OUT/${MID}.hsa.SE.dd.bam  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").sambamba_markdup.log

