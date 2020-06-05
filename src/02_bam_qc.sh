#!/bin/bash

fastq_dir=$1
EXTENSION=$2
BAMDIR=$3
ncells=$( ls ${fastq_dir}/*$EXTENSION | wc -l )
nn=$(( $ncells - 1 ))

bioID=$( basename `pwd` | sed -e 's/_.*//' )

bsub -J fst_${bioID} -n 4 -M 2 -W 359 ./src/qc_scripts/02_flagstat.sh $fastq_dir $EXTENSION $BAMDIR
bsub -J fst_${bioID}[1-$nn] -n 4 -M 2 -W 359 ./src/qc_scripts/02_flagstat.sh $fastq_dir $EXTENSION $BAMDIR

bsub -J idx_${bioID} -n 1 -M 2 -W 359 ./src/qc_scripts/02_idxstats.sh $fastq_dir $EXTENSION $BAMDIR
bsub -J idx_${bioID}[1-$nn] -n 1 -M 2 -W 359 ./src/qc_scripts/02_idxstats.sh $fastq_dir $EXTENSION $BAMDIR

bsub -J sts_${bioID} -n 3 -M 2 -W 359 ./src/qc_scripts/02_stats.sh $fastq_dir $EXTENSION $BAMDIR
bsub -J sts_${bioID}[1-$nn] -n 3 -M 2 -W 359 ./src/qc_scripts/02_stats.sh $fastq_dir $EXTENSION $BAMDIR

bsub -J pcmm_${bioID} -n 1 -M 3 -W 359 ./src/qc_scripts/02_pcmm.sh $fastq_dir $EXTENSION $BAMDIR
bsub -J pcmm_${bioID}[1-$nn] -n 1 -M 3 -W 359 ./src/qc_scripts/02_pcmm.sh $fastq_dir $EXTENSION $BAMDIR

bsub -J qm_${bioID} -n 8 -M 4 -W 359 ./src/qc_scripts/02_qualimap.sh $fastq_dir $EXTENSION $BAMDIR
bsub -J qm_${bioID}[1-$nn] -n 8 -M 4 -W 359 ./src/qc_scripts/02_qualimap.sh $fastq_dir $EXTENSION $BAMDIR

bsub -J precc_${bioID} -n 1 -M 4 -W 359 ./src/qc_scripts/02_pre_ccurve.sh $fastq_dir $EXTENSION $BAMDIR
bsub -J precc_${bioID}[1-$nn] -n 1 -M 4 -W 359 ./src/qc_scripts/02_pre_ccurve.sh $fastq_dir $EXTENSION $BAMDIR

bsub -J ext_${bioID} -n 1 -M 4 -W 359 ./src/qc_scripts/02_pre_extrap.sh $fastq_dir $EXTENSION $BAMDIR
bsub -J ext_${bioID}[1-$nn] -n 1 -M 4 -W 359 ./src/qc_scripts/02_pre_extrap.sh $fastq_dir $EXTENSION $BAMDIR


