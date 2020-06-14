#!/bin/bash

fastq_dir=$1
EXTENSION=$2
ncells=$( ls ${fastq_dir}/*${EXTENSION} | wc -l )
nn=$(( $ncells - 1 ))

bioID=$(basename `pwd` | sed -e 's/_.*//')

bsub -J fqc_${bioID} -n 8 -M 1 -W 89 src/qc_scripts/01_1.fastqc.sh $fastq_dir $EXTENSION
bsub -J fqc_${bioID}[1-$nn] -n 8 -M 1 -W 89 src/qc_scripts/01_1.fastqc.sh $fastq_dir $EXTENSION

bsub -J fp_${bioID} -n 1 -M 6 -W 89 src/qc_scripts/01_1.fastp.sh $fastq_dir $EXTENSION
bsub -J fp_${bioID}[1-$nn] -n 1 -M 6 -W 89 src/qc_scripts/01_1.fastp.sh $fastq_dir $EXTENSION

bsub -J fqs_${bioID} -n 8 -M 8 -W 89 ./src/qc_scripts/01_1.fq_screen.sh $fastq_dir $EXTENSION
bsub -J fqs_${bioID}[1-$nn]%100 -n 8 -M 8 -W 89 ./src/qc_scripts/01_1.fq_screen.sh $fastq_dir $EXTENSION
