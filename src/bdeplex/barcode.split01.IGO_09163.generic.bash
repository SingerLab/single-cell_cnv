#!/bin/bash
set -e -x

./src/bdeplex/barcode.split.gzip.sr03.baslan.py ./btsplit/ ./src/bdeplex/barcode.192.txt <(gunzip -c ./rawdata/JAX_0265/Sample_WD824_P1_IGO_09319_1/WD824_P1_IGO_09319_1_S153_R1_001.fastq.gz ) WD824_P1 > barcode.split.IGO_09319_lane_1.pe03.stdout 2> barcode.split.IGO_09319_lane_1.pe03.stderr &

./src/bdeplex/barcode.split.gzip.sr03.baslan.py ./btsplit/ ./src/bdeplex/barcode.192.txt <(gunzip -c ./rawdata/JAX_0265/Sample_WD824_P6_IGO_09319_2/WD824_P6_IGO_09319_2_S154_R1_001.fastq.gz ) WD824_P6 > barcode.split.IGO_09319_lane_1.pe03.stdout 2> barcode.split.IGO_09319_lane_1.pe03.stderr &
