#!/bin/bash
#BSUB -n 3 -R "rusage[mem=2]" -W 359

## code from mmfansler at Biostars (online at: https://www.biostars.org/p/9864/)
## find $BAM_DIR -name '*.bam' |
##   parallel -j8 -N4095 -m --files samtools merge -u - |
##   parallel --xargs samtools merge -@ 8 merged.bam {}";" rm {}
## BAM_DIR=$1
## find $BAM_DIR -name "*bam" | 
##     parallel --tmpdir ./tmp -j 10 -N284 -m --files samtools merge -u - |
##     parallel --tmpdir ./tmp --xargs samtools merge -@ 10 merged.bam {}";" 


[[ -d tmp/ ]] || mkdir tmp
TMPDIR=tmp/

bam1=$1
bam2=$2
outbam=$3

[[ -f $bam1 ]] || eval 'echo "file $bam1 does not exist" ; exit 1'
[[ -f $bam2 ]] || eval 'echo "file $bam2 does not exist" ; exit 1'

[[ ! -f $outbam ]] || eval 'echo "file $outbam exist" ; exit 1'

BOWTIE_h37=~/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/BowtieIndex/b37.fasta

samtools merge -@ $LSB_MAX_NUM_PROCESSORS --reference $BOWTIE_h37 $outbam $bam1 $bam2
