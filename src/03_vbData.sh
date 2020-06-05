#!/bin/bash
#BSUB -n 3 -M 4 -W 359
echo "aggregating files:"
echo "# number of files in :"
echo "#bin.dir vbData ploidy vbStats" | tr ' ' "\t"
for i in  varbin{5,20,50}k ; do echo ${i}/ $(ls $i/*k50.varbin.data.txt |wc -l ) $( ls $i/*k50.varbin.quantal.stats.txt | wc -l ) $( ls $i/*.varbin.stats.txt | wc -l ) ; done | tr ' ' "\t"

set -e -x -o pipefail -u

MID=${1}
ALIGNER=${2}

[[ ! -z $MID ]] ||  { echo "Sample ID missing!" ; exit 1; }

for i in 50k 20k 5k
do
    Rscript ./src/03_vbData.R  --sample.name=${MID} --input.dir=varbin${i} --output.dir=vbData --bin.size=${i} --aligner=${ALIGNER} &
done

wait ${!}

echo "starting geneCN"
Rscript ./src/04_geneCN.R --sample.name=${MID} --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &
echo "estimating FGA"
Rscript ./src/04_fga.R --sample.name=${MID} --bin.size=5k --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &
echo "simple heatmap"
Rscript ./src/06_vbHeatmap.R --sample.name=${MID} --io.dir=vbData --fig.dir=figures --aligner=${ALIGNER} &

wait ${!}

set +x

echo bioID n.fq.gz n.dd.bam n.bamqc n.fastp n.fastqc n.fastq_screen n.flagstat n.idxstat n.markdup n.aln.metrics n.cycle.metrics n.dist.metrics n.preseq.curve n.presec.extrap n.stats n.5k.quantal n.20k.quantal n.50k.quantal n.mapd.qc n.ploidy n.5k.cells n.20k.cells n.50k.cells | tr ' ' "\t" > processed.files.txt && \
    paste  <( echo $MID )  <(ls bsplit/*gz | wc -l) \
	   <(ls bowtie_out/*dd.bam | wc -l) \
	   <( ls bamqc/ | wc -l ) \
	   <( ls fastp/*ok | wc -l ) \
	   <( ls fastqc/*zip | wc -l ) \
	   <( ls fastq_screen/*html | wc -l ) \
	   <( ls flagstats/ | wc -l ) \
	   <( ls idxstats/ | wc -l ) \
	   <( ls metrics/*markdups | wc -l )  \
	   <( ls metrics/*alignment_summary_metrics | wc -l ) \
	   <( ls metrics/*quality_by_cycle_metrics | wc -l ) \
	   <( ls metrics/*quality_distribution_metrics | wc -l ) \
	   <( ls preseq/*c_curve | wc -l ) \
	   <( ls preseq/*lc_extrap | wc -l ) \
	   <( ls stats/*stats | wc -l) \
	   <( ls varbin5k/*.varbin.data.txt | wc -l ) \
	   <( ls varbin20k/*.varbin.data.txt | wc -l )  \
	   <( ls varbin50k/*.varbin.data.txt | wc -l ) \
	   <( echo $(( $(wc -l vbData/*.50k.k50.varbin.mapd.qc.txt | cut -d ' ' -f 1 ) -1 )) ) \
	   <( echo $(( $(wc -l vbData/*.50k.k50.varbin.quantal.ploidy.txt | cut -d ' ' -f 1 ) -1 )) ) \
	   <( echo $(( $(head -n 1 vbData/*.5k.k50.varbin.data.txt | awk '{print NF}') -3 )) ) \
	   <( echo $(( $(head -n 1 vbData/*.20k.k50.varbin.data.txt | awk '{print NF}') -3 )) ) \
	   <( echo $(( $(head -n 1 vbData/*.50k.k50.varbin.data.txt | awk '{print NF}') -3 )) ) | tr ' ' "\t" >> processed.files.txt 
