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
    echo "bsub -n 8 -M 4 -W 89 ./src/01_bam.dd.sh <genome> <path/to/bam/> <.bam> <bwa_out/> <1.5> <4.8>"
    echo ""
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

echo "LSB_MAX_NUM_PROCESSORS: " $LSB_MAX_NUM_PROCESSORS
echo "MAX_MEM_GB: " $MAX_MEM_GB

GENOME=$1

if [ "$GENOME" = "hsa37" ]; then
    BWA_INDEX=/juno/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
elif [ "$GENOME" = "hsa38" ]; then
    BWA_INDEX=/juno/depot/assemblies/H.sapiens/GRCh38_GDC/index/bwa/0.7.17/GRCh38.d1.vd1.fa
elif [ "$GENOME" = "pdx" ]; then
    BWA_INDEX=/work/singer/genomes/mus_musculus/mm10/Sequence/hg19_mm10_pdx/hg19_mm10.fa
elif [ "$GENOME" = "mmu" ]; then
    BWA_INDEX=/work/singer/genomes/mus_musculus/mm10/Sequence/BWAIndex/GRCm38.p6.fa
fi

## get fastq files
## R1=( $( find -L $1 -name "*R1*" | sort -V ) )
## R2=( $( find -L $1 -name "*R2*" | sort -V ) )
BAM=${2}

EXTENSION=$3 || true
[[ ! -z "$EXTENSION" ]] || EXTENSION=.hsa.PE.md.bam

OUT=$4 || true
[[ ! -z "$OUT" ]] || OUT=bwa_map/
[[ -d $OUT ]] || mkdir $OUT

MID=$( basename ${BAM} $EXTENSION | sed -e 's/_IGO.*//' )
echo $MID

GROUP=$(echo $MID | cut -d "_" -f 1 | sed -E 's/[PRMB]$//')

MIN=${5}
MAX=${6}

## BAM extensions
## PE = Paired-END
## FW = forward pair only
MARKDUP_BAM_EXT=${GENOME}.PE.md.bam
DEDUP_BAM_EXT=${GENOME}.PE.dd.bam
SE_BAM_EXT=${GENOME}.FW.dd.bam

## create log files
[ -d log/${GROUP}/$MID ] || mkdir -p log/${GROUP}/${MID}
[ -d tmp/ ] || mkdir -p tmp/
[ -d metrics/${GROUP} ] || mkdir -p metrics/${GROUP}

## tools
BWA=/work/singer/opt/miniconda3/bin/bwa
SAMTOOLS=/work/singer/opt/miniconda3/bin/samtools
PICARD=/work/singer/opt/miniconda3/bin/picard 
QUALIMAP=/work/singer/opt/miniconda3/bin/qualimap
PRESEQ=/home/gularter/bin/preseq

## if bam file exists: skip to varbin
if [ ! -f  $OUT/${MID}.${SE_BAM_EXT} ] ; then

    echo "## $(date) pipeline start" > ${OUT}/${MID}.ok
    
    ## deduplicate
    $PICARD MarkDuplicates I=${BAM} \
	    REMOVE_DUPLICATES=TRUE \
	    M=metrics/${GROUP}/${MID}.bwa.picard.rmdups \
	    O=$OUT/${MID}.${DEDUP_BAM_EXT} -Xmx$MAX_MEM_GB 2> \
	    log/${GROUP}/${MID}/${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").picard_rmdup.log
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
    ## 1 = genome from arg 1
    ## 2 = bam file
    ## 3 = extension
    ## 4 = min ploidy
    ## 5 = max ploidy
else
    echo "## $(date) $OUT/${MID}.${SE_BAM_EXT} exist, submiiting to varbin + re-QC" > ${OUT}/${MID}.ok
fi
 
bsub -n 1 -M 3 -W 89 ./src/02_varbin_pe.sh $GENOME $OUT/${MID}.${SE_BAM_EXT} ".${SE_BAM_EXT}" ${MIN} ${MAX}

#__end__#


## additional BAM metrics and QC
for i in fastp metrics idxstats stats fastq_screen flagstats preseq bamqc fastqc ; do [ -d $i/${GROUP} ] || mkdir -p ${i}/${GROUP} ; done

## flagstats on all reads

$SAMTOOLS flagstat -@ $LSB_MAX_NUM_PROCESSORS ${BAM} > flagstats/${GROUP}/${MID}.bwa.flagstat  2> log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_flagstat.log
if [ $? -eq 0 ] ; then   echo "## $(date) flagstat complete" >> ${OUT}/${MID}.ok ; else exit 10 ; fi

## all read conting counts
$SAMTOOLS idxstats -@ $LSB_MAX_NUM_PROCESSORS ${BAM} > idxstats/${GROUP}/${MID}.bwa.idxstats  2> log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log
if [ $? -eq 0 ] ; then   echo "## $(date) idxstats complete" >> ${OUT}/${MID}.ok ; else exit 11 ; fi

## unique and FW read conting counts
$SAMTOOLS idxstats -@ $LSB_MAX_NUM_PROCESSORS $OUT/${MID}.${SE_BAM_EXT} > idxstats/${GROUP}/${MID}.FW.bwa.idxstats  2> log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log
if [ $? -eq 0 ] ; then   echo "## $(date) idxstats complete" >> ${OUT}/${MID}.ok ; else exit 11 ; fi

## in all reads of MD file
$SAMTOOLS stats -@ $LSB_MAX_NUM_PROCESSORS ${BAM} > stats/${GROUP}/${MID}.bwa.samtools.stats  2> log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_stats.log
if [ $? -eq 0 ] ; then   echo "## $(date) stats complete" >> ${OUT}/${MID}.ok ; else exit 12 ; fi

## qualimap bamqc all reads
unset DISPLAY
$QUALIMAP bamqc -nt $LSB_MAX_NUM_PROCESSORS -bam ${BAM} -gd HUMAN -outdir bamqc/${GROUP}/${MID}/ -outfile ${MID}.bwa.pdf -outformat "PDF:HTML"  2> log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").qualimap_bamqc.log
if [ $? -eq 0 ] ; then   echo "## $(date) bamqc complete" >> ${OUT}/${MID}.ok ; else exit 13 ; fi

## preseq compexity curve
$PRESEQ c_curve -B ${BAM} > preseq/${GROUP}/${MID}.preseq.c_curve 2> log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.c_curve.log
if [ $? -eq 0 ] ; then   echo "## $(date) preseq c_curve complete" >> ${OUT}/${MID}.ok ; else exit 14 ; fi

## preseq compexity extrapolation
$PRESEQ lc_extrap -B ${BAM} > preseq/${GROUP}/${MID}.preseq.lc_extrap 2> log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.lc_extrap.log
if [ $? -eq 0 ] ; then   echo "## $(date) preseq lc_extrap complete" >> ${OUT}/${MID}.ok ; else exit 15 ; fi

## picard collect multiple metrics
$PICARD CollectMultipleMetrics I=${BAM} O=metrics/${GROUP}/${MID} R=${BWA_INDEX} -Xmx$MAX_MEM_GB  2> log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").picard_metrics.log
if [ $? -eq 0 ] ; then   echo "## $(date) picard cmm complete" >> ${OUT}/${MID}.ok ; else exit 15 ; fi

