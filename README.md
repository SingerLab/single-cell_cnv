# single_cell_cnv
SingerLab pipeline for demultiplexing, processing, and quality control of
fastq files.  Built for use with LSF or direct `bash` command line with LSF environment variables

## dependencies
* Alignment
   - bwa
   - bowtie *version 1*
   - picard
   - samtools
   - sambamba
   
* QC tools
   - fastqc
   - fastq_screen
   - fastp
   - qualimap
   - preseq
   - picard
   - samtools

## Recommended
* QC
   - MultiQC

* Visualization 
   - R/Complex Heatmap
   - IGV

* Analysis
   - R/gac (http://github.com/SingerLab/gac)

## Usage

### steps for single-end 50bp
1. run ./src/bdplex/barcode.bsplit.sh
2. QC fastq files (e.g ./src/01_fastq_qc.sh)
3. run ./src/01_bowtie.sc.map.sh ## automatically launches 02_varbin.sh
4. QC bam files (e.g. ./src/02_bam_qc.sh)
5. run ./src/03_vbData.sh


### steps for paired end 100bp
*Note:* For Paired-end demultiplexing insall wannaAln at https://github.com/aquaskyline/WannaAln

1. run ./src/bdplex/barcode.wsplit.sh
2. QC fastq files
3. run ./src/01_bwa.sc.map.sh
4. QC bam files
5. run 02_varbin.sh
6. run ./src/03_vbData.sh


## Scalability
**Table 1.** Currently processed cells

|bioID    | n.fq.gz| n.dd.bam| n.bamqc| n.fastp| n.fastqc| n.fastq_screen| n.flagstat| n.idxstat| n.markdup| n.aln.metrics| n.cycle.metrics| n.dist.metrics| n.preseq.curve| n.presec.extrap| n.stats| n.5k.quantal| n.20k.quantal| n.50k.quantal| n.mapd.qc| n.ploidy| n.5k.cells| n.20k.cells| n.50k.cells|
|:--------|-------:|--------:|-------:|-------:|--------:|--------------:|----------:|---------:|---------:|-------------:|---------------:|--------------:|--------------:|---------------:|-------:|------------:|-------------:|-------------:|---------:|--------:|----------:|-----------:|-----------:|
|T1       |    1837|     1837|    1837|    1837|     1837|           1549|       1837|      1837|      1837|          1837|            1837|           1837|           1837|            1837|    1837|         1837|          1837|          1837|      1837|     1837|       1837|        1837|        1837|
|T2       |    1968|     1968|    1968|    1968|     1968|            619|       1968|      1968|      1968|          1968|            1968|           1968|           1968|            1968|    1968|         1968|          1968|          1968|      1968|     1968|       1968|        1968|        1968|
|T3       |    1976|     1976|    1976|    1933|     1976|            733|       1976|      1976|      1976|          1976|            1976|           1976|           1976|            1976|    1976|         1976|          1976|          1976|      1976|     1976|       1976|        1976|        1976|
|T4       |     934|      934|     934|     934|      934|            852|        934|       934|       934|           934|             934|            934|            934|             934|     934|          934|           934|           934|       934|      934|        934|         934|         934|
|T5       |    2242|     2242|    2242|    2242|     2242|            966|       2242|      2242|      2242|          2242|            2242|           2242|           2242|            2242|    2242|         2242|          2242|          2242|      2242|     2242|       2242|        2242|        2242|
|T6       |    1056|     1056|    1056|    1056|     1056|            576|       1056|      1056|      1056|          1056|            1056|           1056|           1056|            1056|    1056|         1056|          1056|          1056|      1056|     1056|       1056|        1056|        1056|
|T7       |    1246|     1246|    1246|    1246|     1246|            669|       1246|      1246|      1246|          1246|            1246|           1246|           1246|            1246|    1246|         1246|          1246|          1246|      1246|     1246|       1246|        1246|        1246|
|T8       |    1298|     1298|    1298|    1298|     1298|            132|       1298|      1298|      1298|          1298|            1298|           1298|           1298|            1298|    1298|         1298|          1298|          1298|      1298|     1298|       1298|        1298|        1298|
|T9       |    1604|     1604|    1604|    1604|     1604|           1446|       1604|      1604|      1604|          1604|            1604|           1604|           1604|            1604|    1604|         1604|          1604|          1604|      1604|     1604|       1604|        1604|        1604|
|T10      |    2270|     2270|    2270|    2270|     2270|           2212|       2270|      2270|      2270|          2270|            2270|           2270|           2270|            2270|    2270|         2270|          2270|          2270|      2270|     2270|       2270|        2270|        2270|
|T11      |     290|      290|     290|     290|      290|             27|        290|       290|       290|           286|             286|            286|            290|             290|     290|          290|           290|           290|       290|      290|        290|         290|         290|
|T12      |    1737|     1737|    1737|    1728|     1737|            305|       1737|      1737|      1737|          1737|            1737|           1737|           1737|            1737|    1737|         1737|          1737|          1737|      1737|     1737|       1737|        1737|        1737|
|T13      |     107|      107|     107|     107|      107|              0|        107|       107|       107|           107|             107|            107|            107|             107|     107|          107|           107|           107|       107|      107|        107|         107|         107|
|T14      |      35|       35|      35|      35|       35|              0|         35|        35|        35|            35|              35|             35|             35|              35|      35|           35|            35|            35|        35|       35|         35|          35|          35|
|T15      |     194|      194|     194|     192|      192|              0|        194|       194|       194|           193|             193|            193|            194|             194|     194|          194|           194|           194|       194|      194|        194|         194|         194|
|T16      |    1440|     1440|    1440|    1044|     1056|           1052|       1440|      1440|      1440|          1440|            1440|           1440|           1440|            1440|    1440|         1440|          1440|          1440|      1440|     1440|       1440|        1440|        1440|
