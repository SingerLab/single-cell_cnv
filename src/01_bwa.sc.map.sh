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
    echo " Run times:  for 13M PE reads:"
    echo "  78 min w/ 6 threads and 24G RAM"
    echo "  120 min w/ 4 threads and 24G RAM"
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
MEM_PER_JOB=$(echo $(printf "%d\n" ${LSB_CG_MEMLIMIT} )/1024^3 | bc ) ## $( echo $LSB_SUB_RES_REQ | sed -e 's/.*=//' -e 's/]//' )
MAX_MEM_GB=$(( MEM_PER_JOB*LSB_MAX_NUM_PROCESSORS - 1 ))G

echo "LSB_JOBINDEX: " $LSB_JOBINDEX
echo "LSB_MAX_NUM_PROCESSORS: " $LSB_MAX_NUM_PROCESSORS
echo "MAX_MEM_GB: " $MAX_MEM_GB

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
MARKDUP_BAM_EXT=${GENOME}.PE.md.bam
DEDUP_BAM_EXT=${GENOME}.PE.dd.bam
SE_BAM_EXT=${GENOME}.FW.dd.bam

## create log files
[ -d log/$MID ] || mkdir -p log/$MID
[ -d tmp/ ] || mkdir -p tmp/
[ -d metrics/ ] || mkdir -p metrics/

## if bam file exists: exit
[[ ! -f  $OUT/${MID}.${SE_BAM_EXT} ]] || exit 1

## tools
BWA=/work/singer/opt/miniconda3/bin/bwa
SAMTOOLS=/work/singer/opt/miniconda3/bin/samtools
PICARD=/work/singer/opt/miniconda3/bin/picard 
QUALIMAP=/work/singer/opt/miniconda3/bin/qualimap
PRESEQ=/home/gularter/bin/preseq

echo "## $(date) pipeline start" > ${OUT}/${MID}.ok


## alignment with BWA +- 87 min for 13M reads using 4 cores and 6 of ram
$BWA mem -aM -t $LSB_MAX_NUM_PROCESSORS $BWA_INDEX \
     ${R1[$LSB_JOBINDEX]} ${R2[$LSB_JOBINDEX]} | \
    samtools view -bS - > $OUT/${MID}.bam  2> \
	     log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").bwa.log
if [ $? -eq 0 ] ; then echo "## $(date) alignment complete" >> ${OUT}/${MID}.ok; else exit 2; fi


## collate -- sort by read name
$SAMTOOLS collate -@ $LSB_MAX_NUM_PROCESSORS \
	  -o $OUT/${MID}.cl.bam \
	  $OUT/${MID}.bam
if [ $? -eq 0 ] ; then  echo "## $(date) collate complete" >> ${OUT}/${MID}.ok ; rm ${OUT}/${MID}.bam ; else exit 3 ; fi

$SAMTOOLS fixmate -m -O BAM -@ $LSB_MAX_NUM_PROCESSORS \
	  $OUT/${MID}.cl.bam \
	  $OUT/${MID}.fx.bam
if [ $? -eq 0 ] ; then  echo "## $(date) fixmate complete" >> ${OUT}/${MID}.ok ; rm ${OUT}/${MID}.cl.bam ; else exit 4 ; fi

## sort bam file
$SAMTOOLS sort -@ $LSB_MAX_NUM_PROCESSORS \
 	  -m $MAX_MEM_GB \
 	  -o $OUT/${MID}.sorted.bam \
 	  -O BAM \
 	  $OUT/${MID}.fx.bam  2> \
 	  log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").sort.log
 if [ $? -eq 0 ] ; then  echo "## $(date) sort complete" >> ${OUT}/${MID}.ok ; rm $OUT/${MID}.fx.bam ; else exit 5 ; fi

## remove duplicate reads
$PICARD MarkDuplicates I=$OUT/${MID}.sorted.bam \
	O=$OUT/${MID}.${MARKDUP_BAM_EXT} \
	M=metrics/${MID}.bwa.picard.markdups 2> \
	log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").picard_markdup.log
if [ $? -eq 0 ] ; then  echo "## $(date) mark duplicates complete" >> ${OUT}/${MID}.ok ; rm $OUT/${MID}.sorted.bam ; else exit 6 ; fi

$SAMTOOLS index -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.${MARKDUP_BAM_EXT}

## deduplicate
$PICARD MarkDuplicates I=$OUT/${MID}.${MARKDUP_BAM_EXT} \
	REMOVE_DUPLICATES=TRUE \
	M=metrics/${MID}.bwa.picard.rmdups \
	O=$OUT/${MID}.${DEDUP_BAM_EXT} 2> \
	log/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").picard_rmdup.log
if [ $? -eq 0 ] ; then  echo "## $(date) remove duplicates complete" >> ${OUT}/${MID}.ok ; else exit 7 ; fi

## filter by read quality and unique hits
$SAMTOOLS view -@ $LSB_MAX_NUM_PROCESSORS -q 30 -F 0x800 -o $OUT/${MID}.UQ.bam \
	  $OUT/${MID}.${DEDUP_BAM_EXT}
if [ $? -eq 0 ] ; then   echo "## $(date) filter unique reads complete" >> ${OUT}/${MID}.ok ; rm $OUT/${MID}.${DEDUP_BAM_EXT} ; else exit 8 ; fi

## clip/remove 1 mate of the PE reads (avoids counting twice)
samtools view -@ ${LSB_MAX_NUM_PROCESSORS} -f 0x40 -h -o $OUT/${MID}.${SE_BAM_EXT} \
	 $OUT/${MID}.UQ.bam
if [ $? -eq 0 ] ; then   echo "## $(date) clip PE reads complete" >> ${OUT}/${MID}.ok ; rm $OUT/${MID}.UQ.bam ; else exit 9 ; fi

## indexing
$SAMTOOLS index -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.${SE_BAM_EXT}

## sumbit varbin_pe
bsub -n 1 -M 3 -W 89 ./src/02_varbin_pe.sh $OUT/${MID}.${SE_BAM_EXT} ${SE_BAM_EXT} ${MIN} ${MAX}

#__end__#


## additional BAM metrics and QC
for i in fastp metrics idxstats stats fastq_screen flagstats preseq bamqc fastqc ; do [ -d $i ] || mkdir $i ; done

$SAMTOOLS flagstat -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.${MARKDUP_BAM_EXT} > flagstats/${MID}.bwa.flagstat  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_flagstat.log
if [ $? -eq 0 ] ; then   echo "## $(date) flagstat complete" >> ${OUT}/${MID}.ok ; else exit 10 ; fi

$SAMTOOLS idxstats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.${MARKDUP_BAM_EXT} > idxstats/${MID}.bwa.idxstats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log
if [ $? -eq 0 ] ; then   echo "## $(date) idxstats complete" >> ${OUT}/${MID}.ok ; else exit 11 ; fi

$SAMTOOLS stats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.${MARKDUP_BAM_EXT} > stats/${MID}.bwa.samtools.stats  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_stats.log
if [ $? -eq 0 ] ; then   echo "## $(date) stats complete" >> ${OUT}/${MID}.ok ; else exit 12 ; fi

unset DISPLAY
$QUALIMAP bamqc -nt $LSB_MAX_NUM_PROCESSORS -bam $OUT/${MID}.${MARKDUP_BAM_EXT} -gd HUMAN -outdir bamqc/${MID}/ -outfile ${MID}.bwa.pdf -outformat "PDF:HTML"  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").qualimap_bamqc.log
if [ $? -eq 0 ] ; then   echo "## $(date) bamqc complete" >> ${OUT}/${MID}.ok ; else exit 13 ; fi

$PRESEQ c_curve -B $OUT/${MID}.${MARKDUP_BAM_EXT} > preseq/${MID}.preseq.c_curve 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.c_curve.log
if [ $? -eq 0 ] ; then   echo "## $(date) preseq c_curve complete" >> ${OUT}/${MID}.ok ; else exit 14 ; fi

$PRESEQ lc_extrap -B $OUT/${MID}.${MARKDUP_BAM_EXT} > preseq/${MID}.preseq.lc_extrap 2> log/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.lc_extrap.log
if [ $? -eq 0 ] ; then   echo "## $(date) preseq lc_extrap complete" >> ${OUT}/${MID}.ok ; else exit 15 ; fi

$PICARD CollectMultipleMetrics I=$OUT/${MID}.${MARKDUP_BAM_EXT} O=metrics/$MID R=${BWA_INDEX}  2> log/${MID}/$(date "+%Y%m%d-%H%M%S").picard_metrics.log
if [ $? -eq 0 ] ; then   echo "## $(date) picard cmm complete" >> ${OUT}/${MID}.ok ; else exit 15 ; fi

