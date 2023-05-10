#!/bin/bash
#BSUB -n 8 -R "rusage[mem=5]" -W 359 -R "span[hosts=1]"
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This script to runs the base pipeline to align, FASTQ files to the human"
    echo " GRCh37, mouse mm10, or pdx (GRCh37+mm10) reference genomes.  It then"
    echo " runs mark duplicates and duplicate read removal.  Furthermore, it runs"
    echo "  QC metrics from picard-tools, and samtools on the bam files"
    echo ""
    echo "Usage:"    
    echo "This script expects a genome type <hsa,pdx> as first argument,  directory"
    echo " name with the FASTQ files as the second argument; the FASTQ extension as"
    echo " the third e.g. _001.fastq.gz, .fq.gz; the output directory name as the"
    echo " fourth argument; and the ploidy range as \`min' \`max' on the fifth and"
    echo " sixth arguments, respectively."
    echo ""
    echo "Example:"
    echo "bsub -n 8 -M 4 -W 89 ./src/01_bwa.sc.map.sh <genome> <path/to/fastq/> <.fq.gz> <bwa_out/> <1.5> <4.8>"
    echo "bsub -J 'bm[1-1500]' -n 8 -M 4 -W 89 ./src/01_bwa.sc.map.sh hsa path/to/fastq/ .fq.gz bwa_out/ 1.5 4.8"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 8 -M 4 -W 89 LSB_JOBINDEX=9 ./src/01_bwa.sc.map.sh hsa path/to/fastq/ .fq.gz bwa_out/ 1.5 4.8"
    echo ""
    exit 1;
}
#</usage>

## suggestions from
## https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -e -x -E -u -o pipefail

## trap read debug
MAX_MEM_GB=38G

echo "LSB_JOBINDEX: " $LSB_JOBINDEX
echo "LSB_MAX_NUM_PROCESSORS: " $LSB_MAX_NUM_PROCESSORS

GENOME=$1

if [ "$GENOME" = "hsa" ]; then
    BWA_INDEX=/juno/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
elif [ "$GENOME" = "pdx" ]; then
    BWA_INDEX=/work/singer/genomes/mus_musculus/mm10/Sequence/hg19_mm10_pdx/hg19_mm10.fa
elif [ "$GENOME" = "mmu" ]; then
    BWA_INDEX=/work/singer/genomes/mus_musculus/mm10/Sequence/BWAIndex/GRCm38.p6.fa
fi

## get fastq files
## R1=( $( find -L $1 -name "*R1*" | sort -V ) )
## R2=( $( find -L $1 -name "*R2*" | sort -V ) )
R1=( $( ls ${2}*R1.fq.gz | sort -V ) )
R2=( $( ls ${2}*R2.fq.gz | sort -V ) )
echo ${R1[@]}
echo ${R2[@]}

EXTENSION=$3 || true
[[ ! -z "$EXTENSION" ]] || EXTENSION=.R1.fq.gz

OUT=$4 || true
[[ ! -z "$OUT" ]] || OUT=bwa_map/
[[ -d $OUT ]] || mkdir $OUT

MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//' )

echo $MID

## BAM extensions
## PE = Paired-END
## ${GENOME}.PE.dd.bam -- change to choose genomes
DEDUP_BAM_EXT=${GENOME}.PE.dd.bam
MARKDUP_BAM_EXT=${GENOME}.PE.md.bam

## create log files
[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/
[ -d metrics/ ] || mkdir -p metrics/

## if bam file exists: exit
[[ ! -f  $OUT/${MID}.${MARKDUP_BAM_EXT} ]] || exit 1

## tools
BWA=/work/singer/opt/miniconda3/bin/bwa
SAMTOOLS=/work/singer/opt/miniconda3/bin/samtools
SAMBAMBA=/home/gularter/bin/sambamba
PICARD=/work/singer/opt/miniconda3/bin/picard 
QUALIMAP=/work/singer/opt/miniconda3/bin/qualimap
PRESEQ=/home/gularter/bin/preseq

## alignment with BWA
$BWA mem -aM -t $LSB_MAX_NUM_PROCESSORS $BWA_INDEX \
     ${R1[$LSB_JOBINDEX]} ${R2[$LSB_JOBINDEX]} | \
    samtools view -bS - > $OUT/${MID}.bam  2> \
	     log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").bwa.log
if [ $? -eq 0 ] ; then echo "alignment succesful" ; else exit 1; fi

## sort bam file
$SAMBAMBA sort -t $LSB_MAX_NUM_PROCESSORS \
	  -m $MAX_MEM_GB \
	  --tmpdir=tmp/ \
	  $OUT/${MID}.bam  2> \
	  log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").sort.log
if [ $? -eq 0 ] ; then touch ${OUT}/${MID}.bam.ok ; else exit 2 ; fi

## remove duplicate reads
$PICARD MarkDuplicates I=$OUT/${MID}.sorted.bam \
	O=$OUT/${MID}.${MARKDUP_BAM_EXT} \
	M=metrics/${MID}.bwa.picard.markdups  2> \
	log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").picard_markdup.log
if [ $? -eq 0 ] ; then rm $OUT/${MID}.sorted.bam{,.bai} ${OUT}/${MID}.bam{,.ok} ; else exit 6 ; fi

## indexing
$SAMTOOLS index -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.${MARKDUP_BAM_EXT}


#__end__#


## additional metrics and QC
#% for i in fastp metrics idxstats stats fastq_screen flagstats preseq bamqc fastqc ; do [ -d $i ] || mkdir $i ; done
#% 
#% $SAMTOOLS flagstat -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.md.bam > flagstats/${MID}.bowtie.flagstat  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_flagstat.log
#% if [ $? -eq 0 ] ; then echo "samtools flagstat succesfull" ; else exit 7 ; fi
#% 
#% $SAMTOOLS idxstats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.md.bam > idxstats/${MID}.bowtie.idxstats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log
#% if [ $? -eq 0 ] ; then echo "samtools idxstats succesfull" ; else exit 8 ; fi
#% 
#% $SAMTOOLS stats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.md.bam > stats/${MID}.bowtie.samtools.stats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_stats.log
#% if [ $? -eq 0 ] ; then echo "samtools stats succesfull" ; else exit 9 ; fi
#% 
#% unset DISPLAY
#% $QUALIMAP bamqc -nt $LSB_MAX_NUM_PROCESSORS -bam $OUT/${MID}.md.bam -gd HUMAN -outdir bamqc/${MID}/ -outfile ${MID}.bowtie.pdf -outformat "PDF:HTML"  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").qualimap_bamqc.log
#% if [ $? -eq 0 ] ; then echo "qualimap succesfull" ; else exit 10 ; fi
#% 
#% $PRESEQ c_curve -B $OUT/${MID}.md.bam > preseq/${MID}.preseq.c_curve 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.c_curve.log
#% if [ $? -eq 0 ] ; then echo "preseq c_curve succesfull" ; else exit 11 ; fi
#% 
#% $PRESEQ lc_extrap -B $OUT/${MID}.md.bam > preseq/${MID}.preseq.lc_extrap 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.lc_extrap.log
#% if [ $? -eq 0 ] ; then echo "preseq c_curve succesfull" ; else exit 12 ; fi
#% 
#% ## placed here as it often fails
#% $PICARD CollectMultipleMetrics I=$OUT/${MID}.md.bam O=metrics/$MID R=${BOWTIE_h37}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").picard_metrics.log
#% if [ $? -eq 0 ] ; then echo "picard metrics succesfull" ; else exit 13 ; fi
#% 
