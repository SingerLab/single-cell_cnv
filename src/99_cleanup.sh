#!/bin/bash

[[ $# -gt 0 ]] || {
    echo "Usage:"
    echo "requires project input directory e.g. "
    echo "./src/99_clenaup.sh WD5816_igo_09357/"
    exit 1;
}

echo "WARNING! You are deliting tmp files, qc files, and data that is easily generated, or stored elsewhere"

INDIR=$1

echo "cleanup start: " `date` >> ${INDIR}/99_cleanup.log

find ${INDIR}/seqdata/ -name "*Sample*IGO*fastq.gz" -exec rm {} \;
echo "fastq files have been removed by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/bowtie_out/ -name "*md.bam" -exec rm {} \;
echo "*md.bam files have been removed by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/bowtie_out/ -name "*md.bam.bai" -exec rm {} \;
echo "*md.bam.bai files have been removed by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/bowtie_out/ -name "*ok" -exec rm {} \;
echo "bowtie/*.ok files have been removed by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

## find ${INDIR}/bsplit/low_counts/ -exec rm {} \;
find ${INDIR}/bsplit/unexpected/ -exec rm {} \;
echo "low_counts/ and unexpected/ have been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/bamqc/ -exec rm {} \;
find ${INDIR}/bamqc/ -exec rmdir {} \;
find ${INDIR}/bamqc/ -exec rmdir {} \;
echo "bamqc/ has been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/fastp/ -exec rm {} \;
echo "fastp/ has been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/fastqc/ -exec rm {} \;
echo "fastqc/ has been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/fastq_screen/ -name "*fq.gz_temp_subset.fastq" -exec rm {} \;
echo "*fq.gz_temp_subset.fastq have been removed by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/fastq_screen/ -name "*aligner_standard_error*" -exec rm {} \; 
echo "*aligner_standard_error* have been removed by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/fastq_screen/ -name "*html" -exec rm {} \;
echo "fastq_screen/ has been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/log/ -exec rm {} \;
echo "log/ has been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/metrics/ -name "*pdf" -exec rm {} \;
echo "metrics/*pdf have been removed by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/preseq/ -exec rm {} \;
echo "preseq/ has been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/stats/ -exec rm {} \;
echo "stats/ has been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/tmp/ -exec rm {} \;
find ${INDIR}/tmp/ -exec rmdir {} \;
echo "tmp/ has been cleared by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

find ${INDIR}/varbin*/ -name '*dd.*.varbin.out.txt' -exec rm {} \;
find ${INDIR}/varbin*/ -name '*dd.*.varbin.stats.txt' -exec rm {} \;
find ${INDIR}/varbin*/ -name '*dd.grch37.*.k50*varbin.short.txt' -exec rm {} \;
find ${INDIR}/varbin*/ -name '*dd.grch37.*.k50*varbin.short.cbs.seg' -exec rm {} \;
find ${INDIR}/varbin*/ -name '*dd.grch37.*.k50*varbin.quantal.stats.txt' -exec rm {} \;
find ${INDIR}/varbin*/ -name '*dd.grch37.*.k50.varbin.data.txt' -exec rm {} \;
find ${INDIR}/varbin*/ -name '*dd.5k.wg.png' -exec rm {} \;
find ${INDIR}/varbin*/ -name '*dd.5k.wg.quantal.png' -exec rm {} \;
find ${INDIR}/varbin*/ -name '*dd.5k.wg.nobad.png' -exec rm {} \;
echo "varbin5k/ varbin20k/ and varbin50/ have been cleaned by 99_cleanup.sh" >> ${INDIR}/99_cleanup.log

echo "cleanup end: " `date` >> ${INDIR}/99_cleanup.log
