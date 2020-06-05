#!/bin/bash
#BSUB -n 10 -R "rusage[mem=2]" -W 359

## code from mmfansler at Biostars (online at: https://www.biostars.org/p/9864/)
## find $BAM_DIR -name '*.bam' |
##   parallel -j8 -N4095 -m --files samtools merge -u - |
##   parallel --xargs samtools merge -@8 merged.bam {}";" rm {}

[[ -d tmp/ ]] || mkdir tmp
TMPDIR=tmp/

BAM_DIR=$1
find $BAM_DIR -name "*bam" | 
    parallel --tmpdir ./tmp -j 10 -N284 -m --files samtools merge -u - |
    parallel --tmpdir ./tmp --xargs samtools merge -@ 10 merged.bam {}";" 

