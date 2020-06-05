#!/bin/bash
set -x

BWGA192_IDX=( $( cut -f 2 $1 ) )
BWGA192_SAMPLE=( $( cut -f 3 $1 ) )

N_SAMPLES=${#BWGA192_SAMPLE[@]}

echo ${BWGA192_IDX[@]}
echo ${BWGA192_SAMPLE[@]}

SUB=$2

for i in {1..61}
do    
    echo "renaming ${SUB}${BWGA192_IDX[$i]}.fastq.gz to ${BWGA192_SAMPLE[$i]}.${BWGA192_IDX[$i]}.fastq.gz"
    mv ${SUB}${BWGA192_IDX[$i]}.fq.gz ${BWGA192_SAMPLE[$i]}.${BWGA192_IDX[$i]}.fq.gz
done

	 
## for i in {1..61}
## do    
##     echo "renaming ${SUB}${BWGA192_IDX[$i]}.R1.fastq.gz to ${BWGA192_SAMPLE[$i]}.${BWGA192_IDX[$i]}.R1.fastq.gz"
##     mv ${SUB}${BWGA192_IDX[$i]}.fastq.gz ${BWGA192_SAMPLE[$i]}.${BWGA192_IDX[$i]}.fastq.gz
##     echo "renaming ${SUB}${BWGA192_IDX[$i]}.R2.fastq.gz to ${BWGA192_SAMPLE[$i]}.${BWGA192_IDX[$i]}.R2.fastq.gz"
##     mv ${SUB}${BWGA192_IDX[$i]}.R2.fastq.gz ${BWGA192_SAMPLE[$i]}.${BWGA192_IDX[$i]}.R2.fastq.gz
## done

