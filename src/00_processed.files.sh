#!/bin/bash

echo bioID n.fq.gz n.dd.bam n.md.bam n.bamqc n.fastp n.fastqc n.fastq_screen n.flagstat n.idxstat n.markdup n.aln.metrics n.cycle.metrics n.dist.metrics n.preseq.curve n.presec.extrap n.stats n.5k.quantal n.20k.quantal n.50k.quantal n.mapd.qc n.ploidy n.5k.cells n.20k.cells n.50k.cells | tr ' ' "\t"

paste  <( echo $1 )  <(ls $2/*gz | wc -l) \
       <(ls $3/*dd.bam | wc -l) \
       <(ls $3/*md.bam | wc -l) \
       <( ls bamqc/ | wc -l ) \
       <( ls fastp/*html | wc -l ) \
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
       <( ls varbin5k/*k50.varbin.data.txt | wc -l ) \
       <( ls varbin20k/*k50.varbin.data.txt | wc -l )  \
       <( ls varbin50k/*k50.varbin.data.txt | wc -l ) \
       <( echo $(( $(wc -l vbData/*.50k.k50.varbin.mapd.qc.txt | cut -d ' ' -f 1 ) -1 )) ) \
       <( echo $(( $(wc -l vbData/*.50k.k50.varbin.quantal.ploidy.txt | cut -d ' ' -f 1 ) -1 )) ) \
       <( echo $(( $(zcat vbData/*.5k.k50.varbin.data.txt.gz | head -n 1 | awk '{print NF}') -3 )) ) \
       <( echo $(( $(zcat vbData/*.20k.k50.varbin.data.txt.gz | head -n 1 | awk '{print NF}') -3 )) ) \
       <( echo $(( $(zcat vbData/*.50k.k50.varbin.data.txt.gz | head -n 1 | awk '{print NF}') -3 )) ) | tr ' ' "\t"
