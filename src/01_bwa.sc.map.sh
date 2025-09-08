#!/bin/bash
#BSUB -n 8 -R "rusage[mem=5]" -W 359 -R "span[hosts=1]"

#<usage>
# Usage function
usage() {
    echo "Description:"
    echo ""
    echo "This script to runs the base pipeline to align, FASTQ files to human,"
    echo " mouse, or pdx (GRCh37+mm10) reference genomes.  It then"
    echo " runs mark duplicates and duplicate read removal.  Furthermore, it runs"
    echo "  QC metrics from picard-tools, and samtools on the bam files"
    echo ""
    echo "Usage: $0 -g GENOME [options] /path/to/fastq/"
    echo ""
    echo "Options:"
    echo "  -g,  --genome         Path to reference genome (required)"
    echo "                        accepted genomes:  hsa37 hsa38 mm10 mm38 mm39 pdx37"
    echo "  -e,  --extension      FASTQ file extension (default: _R1.fastq.gz)"
    echo "  -o,  --output_dir     Directory to save outputs (default: ./bwa_out/)"
    echo "  -m,  --min_ploidy     Minimum ploidy (default: 1.5)"
    echo "  -x,  --max_ploidy     Maximum ploidy (default: 5.5)"
    echo "  -h,  --help           Display help message"
    echo ""
    echo " Run times:  for 13M PE reads:"
    echo "  78 min w/ 6 threads and 4G RAM/thread"
    echo "  120 min w/ 4 threads and 6G RAM/thread"
    echo ""
    echo "Example:"
    echo "bsub -J 'bm[1-1500]%200' -n 8 -M 4 -W 89 ./src/01_bwa.sc.map.sh -g hsa38 -e _R1.fq.gz -o bwa_out/ -m 1.5 -x 4.8 path/to/fastq"
    echo ""
    echo "Alternatively it accepts LSB_JOBINDEX as an environment variable to run"
    echo " specific files"
    echo "bsub -n 8 -M 4 -W 89 LSB_JOBINDEX=9 ./src/01_bwa.sc.map.sh -g hsa38 -e _R1.fq.gz -o bwa_out/ -m 1.5 -x 4.8 path/to/fastq "
    echo ""
	
  exit 1
}
#</usage>


# Default values
EXTENSION="_R1.fastq.gz"
OUTPUT_DIR="./bwa_out/"
MIN_PLOIDY=1.5
MAX_PLOIDY=5.5

# Parse options
ARGS=$(getopt -o g:e:o:m:x: -l genome:,extension:,output_dir:,min_ploidy:,max_ploidy:,help -n "$0" -- "$@")
if [ $? -ne 0 ]; then usage; fi

# Reorganize the arguments as parsed by getopt
eval set -- "$ARGS"

while true; do
  case "$1" in
    -g|--genome) GENOME="$2"; shift 2;;
    -e|--extension) EXTENSION="$2"; shift 2;;
    -o|--output_dir) OUTPUT_DIR="$2"; shift 2;;
    -m|--min_ploidy) MIN_PLOIDY="$2"; shift 2;;
    -x|--max_ploidy) MAX_PLOIDY="$2"; shift 2;;
    -h|--help) usage;;
    --) shift; break;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

# Expect FASTQ_DIR as a positional argument
if [[ -z "$1" ]]; then
  echo "Error: FASTQ_DIR (positional) is required"
  usage
fi
FASTQ_DIR="$1"
shift

# Check required inputs
if [[ -z "$GENOME" ]]; then
  echo "Error: --genome is required"
  usage
fi

# Echo settings for debug
echo "Genome: $GENOME"
echo "Extension: $EXTENSION"
echo "Output Dir: $OUTPUT_DIR"
echo "Min Ploidy: $MIN_PLOIDY"
echo "Max Ploidy: $MAX_PLOIDY"
echo "FASTQ_DIR" $FASTQ_DIR

# === Processing ===
## suggestions from
## https://vaneyckt.io/posts/saferd_bash_scripts_with_set_euxo_pipefail/
set -e -x -E -u -o pipefail

## trap read debug
MEM_PER_JOB=$(echo $(printf "%d\n" ${LSB_CG_MEMLIMIT} )/1024^3 | bc )
MAX_MEM_GB=$(( MEM_PER_JOB*LSB_MAX_NUM_PROCESSORS - 1 ))G

echo "LSB_JOBINDEX: " $LSB_JOBINDEX
echo "LSB_MAX_NUM_PROCESSORS: " $LSB_MAX_NUM_PROCESSORS
echo "MAX_MEM_GB: " $MAX_MEM_GB


## conditional genome references
if [ "$GENOME" = "hsa37" ]; then
    BWA_INDEX=/juno/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
elif [ "$GENOME" = "hsa38" ]; then
    BWA_INDEX=/work/singer/genomes/homo_sapiens/Ensembl/GRCh38.p5/Sequence/BWAIndex/GRCh38.d1.vd1.fa.gz
elif [ "$GENOME" = "pdx37" ]; then
    BWA_INDEX=/work/singer/genomes/mus_musculus/mm10/Sequence/hg19_mm10_pdx/hg19_mm10.fa
elif [ "$GENOME" = "mm10" ]; then
    BWA_INDEX=/work/singer/genomes/mus_musculus/mm10/Sequence/BWAIndex/GRCm38.p6.fa
elif [ "$GENOME" = "mm38" ]; then
    BWA_INDEX=/work/singer/genomes/mus_musculus/mm10/Sequence/BWAIndex/GRCm38.p6.fa
elif [ "$GENOME" = "mm39" ]; then
    BWA_INDEX=/work/singer/genomes/mus_musculus/mm39/Sequence/BWAIndex/mm39.fa    
fi

## get fastq file array
R1=( $( ls ${FASTQ_DIR}*${EXTENSION} | sort -V ) )
## R2=( $( ls ${FASTQ_DIR}*R2.fq.gz | sort -V ) )
R2=( $( for i in ${R1[@]} ; do echo $i | sed 's/_R1/_R2/' ; done ) ) 

## check if it exists or create
[[ -d $OUTPUT_DIR ]] || mkdir $OUTPUT_DIR

## extract base name of sample
MID=$( basename ${R1[$LSB_JOBINDEX]} $EXTENSION | sed -e 's/_IGO.*//' )

## extract group / bioID
GROUP=$(basename $FASTQ_DIR | cut -d "_" -f 1 | sed -E 's/[PRMB]$//')

## print group and sample ID
echo $GROUP $MID

## BAM extensions
## PE = Paired-END
## ${GENOME}.PE.dd.bam -- change to choose genomes
MARKDUP_BAM_EXT=${GENOME}.PE.md.bam
DEDUP_BAM_EXT=${GENOME}.PE.dd.bam
SE_BAM_EXT=${GENOME}.FW.dd.bam

## create log directories
[ -d log/${GROUP}/${MID} ] || mkdir -p log/${GROUP}/${MID}/
[ -d tmp/${GROUP} ] || mkdir -p tmp/${GROUP}/

## create qc_stats directory
QC_STATS=qc_stats/

## qc directory setup
for i in metrics ; do 
    [ -d ${QC_STATS}/${GROUP}/${i}/ ] || mkdir -p ${QC_STATS}/${GROUP}/${i}/
done


## no longer needed as alignment is checkpointed
## if bam file exists: exit
## [[ ! -f  ${OUTPUT_DIR}/${MID}.${SE_BAM_EXT} ]] || exit 1

## tools
BWA=/work/singer/opt/miniconda3/bin/bwa
SAMTOOLS=/work/singer/opt/miniconda3/bin/samtools
PICARD=/work/singer/opt/miniconda3/bin/picard 
QUALIMAP=/work/singer/opt/miniconda3/bin/qualimap
PRESEQ=/home/gularter/bin/preseq
FASTQC=/work/singer/opt/miniconda3/bin/fastqc

## checkpoint
OK_FILE=${OUTPUT_DIR}/${MID}.${GENOME}.ok

# Function to update checkpoint in OK file
update_checkpoint() {
    echo "## $(date) $1 complete" >> "$OK_FILE"
}

# Function to check if step has been completed
check_checkpoint() {
    grep -q "$1 complete" "$OK_FILE" 2> /dev/null
}

# Pipeline start
if [ ! -f $OK_FILE ] ; then
    echo "## $(date) pipeline start" > "$OK_FILE"
fi


## alignment with BWA +- 87 min for 13M reads using 4 cores and 6 of ram
if ! check_checkpoint "alignment"; then
    $BWA mem -aM -t $LSB_MAX_NUM_PROCESSORS $BWA_INDEX \
	 ${R1[$LSB_JOBINDEX]} ${R2[$LSB_JOBINDEX]} 2> \
	 log/${GROUP}/${MID}/${MID}_${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").bwa.log | \
	samtools view -bS - > ${OUTPUT_DIR}/${MID}.bam

    if [ $? -eq 0 ] ; then
	update_checkpoint "alignment"
    else
	exit 2
    fi
fi

	
## collate -- sort by read name
if ! check_checkpoint "collate"; then
    $SAMTOOLS collate -@ $LSB_MAX_NUM_PROCESSORS \
	      -o ${OUTPUT_DIR}/${MID}.cl.bam \
	      ${OUTPUT_DIR}/${MID}.bam
    if [ $? -eq 0 ] ; then
	update_checkpoint "collate"
	[ -f ${OUTPUT_DIR}/${MID}.bam ] && \
	    rm ${OUTPUT_DIR}/${MID}.bam
    else
	exit 3
    fi
fi


## fixmate
if ! check_checkpoint "fixmate"; then
    $SAMTOOLS fixmate -m -O BAM -@ $LSB_MAX_NUM_PROCESSORS \
	      ${OUTPUT_DIR}/${MID}.cl.bam \
	      ${OUTPUT_DIR}/${MID}.fx.bam
    if [ $? -eq 0 ]; then
	update_checkpoint "fixmate"
	[ -f ${OUTPUT_DIR}/${MID}.cl.bam ] && \
	    rm ${OUTPUT_DIR}/${MID}.cl.bam
    else
	exit 4
    fi
fi


## sort bam file
if ! check_checkpoint "sort"; then
    $SAMTOOLS sort -@ $LSB_MAX_NUM_PROCESSORS \
 	      -m ${MEM_PER_JOB}g \
 	      -o ${OUTPUT_DIR}/${MID}.sorted.bam \
 	      -O BAM \
 	      ${OUTPUT_DIR}/${MID}.fx.bam  2> \
 	      log/${GROUP}/${MID}/${MID}_${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").sort.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "sort"
	[ -f ${OUTPUT_DIR}/${MID}.fx.bam ] && \
	    rm ${OUTPUT_DIR}/${MID}.fx.bam
    else
	exit 5
    fi
fi


## mark duplicate reads
if ! check_checkpoint "mark duplicates"; then
    $PICARD MarkDuplicates I=${OUTPUT_DIR}/${MID}.sorted.bam \
	    O=${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} \
	    M=${QC_STATS}/${GROUP}/metrics/${MID}.bwa.picard.markdups 2> \
	    log/${GROUP}/${MID}/${MID}_${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").picard_markdup.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "mark duplicates"
	[ -f ${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} ] && \
	    rm ${OUTPUT_DIR}/${MID}.sorted.bam
    else
	exit 6
    fi
fi


if ! check_checkpoint "index MD bam"; then
    $SAMTOOLS index -@ $LSB_MAX_NUM_PROCESSORS \
	      ${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT}
        if [ $? -eq 0 ] ; then
	update_checkpoint "index MD bam"
    else
	exit 6
    fi
fi


## remove duplicate reads 
if ! check_checkpoint "deduplicate"; then
    $PICARD MarkDuplicates I=${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} \
	REMOVE_DUPLICATES=TRUE \
	M=${QC_STATS}/${GROUP}/metrics/${MID}.bwa.picard.rmdups \
	O=${OUTPUT_DIR}/${MID}.${DEDUP_BAM_EXT} 2> \
	log/${GROUP}/${MID}/${MID}_${LSB_JOBID}_$(date "+%Y%m%d-%H%M%S").picard_rmdup.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "deduplicate"
    else
	exit 7
    fi
fi


## filter by read quality and unique hits
if ! check_checkpoint "filter unique reads"; then
    $SAMTOOLS view -@ $LSB_MAX_NUM_PROCESSORS -q 30 -F 0x800 \
	      -o ${OUTPUT_DIR}/${MID}.UQ.bam \
	      ${OUTPUT_DIR}/${MID}.${DEDUP_BAM_EXT}
    if [ $? -eq 0 ] ; then
	update_checkpoint "filter unique reads"
	[ -f ${OUTPUT_DIR}/${MID}.${DEDUP_BAM_EXT} ]  && \
	    rm ${OUTPUT_DIR}/${MID}.${DEDUP_BAM_EXT}
    else
	exit 8
    fi
fi


## clip/remove 1 mate of the PE reads (avoids counting twice)
if ! check_checkpoint "clip PE reads"; then 
    samtools view -@ ${LSB_MAX_NUM_PROCESSORS} -f 0x40 -h \
	     -o ${OUTPUT_DIR}/${MID}.${SE_BAM_EXT} \
	     ${OUTPUT_DIR}/${MID}.UQ.bam
    if [ $? -eq 0 ] ; then
	update_checkpoint "clip PE reads"
	[ -f ${OUTPUT_DIR}/${MID}.UQ.bam ] && 
	    rm ${OUTPUT_DIR}/${MID}.UQ.bam
    else
	exit 9
    fi
fi

    
## indexing
if ! check_checkpoint "index SE bam"; then
    $SAMTOOLS index -@ $LSB_MAX_NUM_PROCESSORS \
	      ${OUTPUT_DIR}/${MID}.${SE_BAM_EXT}
    if [ $? -eq 0 ] ; then
	update_checkpoint "index SE bam"
    else
	exit 10
    fi
fi

## sumbit varbin_pe
if ! check_checkpoint "varbin_pe"; then
    bsub -n 1 -M 3 -W 89 ./src/02_varbin_pe.sh $GENOME \
	 ${OUTPUT_DIR}/${MID}.${SE_BAM_EXT} \
	 ${SE_BAM_EXT} \
	 ${MIN_PLOIDY} ${MAX_PLOIDY}
    if [ $? -eq 0 ] ; then
	update_checkpoint "varbin_pe"
    fi
fi


## __end_of_pipeline__

## __estimate_qc_metrics_start__

## qc directory setup
for i in metrics idxstats stats fastq_screen flagstats preseq bamqc fastqc ; do
    [ -d ${QC_STATS}/${GROUP}/${i}/ ] || mkdir -p ${QC_STATS}/${GROUP}/${i}/
done

## samtools flagstat
if ! check_checkpoint "flagstat"; then
    $SAMTOOLS flagstat -@ $LSB_MAX_NUM_PROCESSORS ${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} > \
	      ${QC_STATS}/${GROUP}/flagstats/${MID}.${GENOME}.bwa.flagstat  2> \
	      log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_flagstat.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "flagstat"
    else
	exit 11
    fi
fi


## samtools index stats
if ! check_checkpoint "idxstats"; then
    $SAMTOOLS idxstats -@ $LSB_MAX_NUM_PROCESSORS ${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} > \
	      ${QC_STATS}/${GROUP}/idxstats/${MID}.${GENOME}.bwa.idxstats  2> \
	      log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_idxstats.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "idxstats"
    else
	exit 12
    fi
fi



## samtools stats
if ! check_checkpoint "samtools stats"; then
$SAMTOOLS stats -@ $LSB_MAX_NUM_PROCESSORS ${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} > \
	  ${QC_STATS}/${GROUP}/stats/${MID}.${GENOME}.bwa.samtools.stats  2> \
	  log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").samtools_stats.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "samtools stats"
    else
	exit 13
    fi
fi



## bamqc qualimap
if ! check_checkpoint "bamqc"; then
    unset DISPLAY
    $QUALIMAP bamqc -nt $LSB_MAX_NUM_PROCESSORS \
	      -bam ${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} \
	      -gd HUMAN -outdir ${QC_STATS}/${GROUP}/bamqc/${MID}/ \
	      -outfile ${MID}.${GENOME}.bwa.pdf -outformat "PDF:HTML"  2> \
	      log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").qualimap_bamqc.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "bamqc"
    else
	exit 14
    fi
fi



## preseq complexity curve
if ! check_checkpoint "preseq c_curve" ; then
    $PRESEQ c_curve -B ${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} > \
	${QC_STATS}/${GROUP}/preseq/${MID}.${GENOME}.preseq.c_curve 2> \
	log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.c_curve.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "preseq c_curve"
    else
	exit 15
    fi
fi


## fastqc
if ! check_checkpoint "preseq lc_extrap" ; then
    $PRESEQ lc_extrap -B ${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} > \
	    ${QC_STATS}/${GROUP}/preseq/${MID}.${GENOME}.preseq.lc_extrap 2> \
	    log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").preseq.lc_extrap.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "preseq lc_extrap"
    else
	exit 16
    fi
fi


## picard collect multiple metrics
if ! check_checkpoint "picard cmm" ; then
    $PICARD CollectMultipleMetrics I=${OUTPUT_DIR}/${MID}.${MARKDUP_BAM_EXT} \
	    O=${QC_STATS}/${GROUP}/metrics/${MID}.${GENOME} R=${BWA_INDEX}  2> \
	    log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").picard_metrics.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "picard cmm"
    else
	exit 17
    fi
fi

## fastqc
if ! check_checkpoint "fastqc" ; then
    $FASTQC -t 2 --outdir ${QC_STATS}/${GROUP}/fastqc/ \
	    ${R1[$LSB_JOBINDEX]} ${R2[$LSB_JOBINDEX]} 2> \
	    log/${GROUP}/${MID}/$(date "+%Y%m%d-%H%M%S").fastqc.log
    if [ $? -eq 0 ] ; then
	update_checkpoint "fastqc"
    else
	exit 18
    fi
fi

## submit fastq_screen
if ! check_checkpoint "fastq_screen" ; then
    bsub -n 8 -M 6 -W 24:00 src/qc_scripts/01_1.fq_screen.sh \
	 $FASTQ_DIR/${MID} ${EXTENSION} ${QC_STATS}/${GROUP}/fastq_screen/
    
    if [ $? -eq 0 ] ; then
	update_checkpoint "fastq_screen"
    else
	exit 19
    fi
fi

## __EOF__

