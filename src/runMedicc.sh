#!/bin/bash

set -e -x

SAMPLE=$1
DESC=$2
OUTDIR=${1:-$3}

for i in $( ls medicc/${SAMPLE}.chr*fa | cut -d '.' -f 2 | sort -V | uniq ) ; do echo $i ${SAMPLE}.${i}.major.cnv.fa ${SAMPLE}.${i}.minor.cnv.fa ; done | tr ' ' "\t" > ${DESC}

## -v = verbose ; -d name of diploid cell ; -c asume molecular clock (Fitch Margoliash)
~/src/medicc/medicc.py -v -d "diploid" ${DESC} ${OUTDIR}
