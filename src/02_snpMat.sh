#!/bin/bash
#BSUB -n 1 -M 4 -W 359

set -e -x

BWA_h37=/ifs/depot/pi/resources/genomes/GRCh37/bwa_fasta/b37.fasta
GFF_h37=/ifs/work/socci/Work/Users/SingerS/GularteR/genomes/homo_sapiens/Ensembl/GRCh37.p13/Annotation/Genes/gencode.v19.annotation.gtf
DBSNP=/ifs/work/socci/Work/Users/SingerS/GularteR/genomes/homo_sapiens/Ensembl/GRCh37.p13/Annotation/Variation/1000genomes/1000GENOMES-phase_3.sorted.vcf.gz

source /home/gularter/opt/miniconda3/bin/activate dna-cnv
PATH=~/src/samtools/:$PATH

OUTDIR=facets/snp_pileup
[[ -d $OUTDIR ]] || mkdir -p $OUTDIR

TUMOR=$1
NORMAL=$2

TM=$( basename $TUMOR .dd.bam )
NN=$( basename $NORMAL .dd.bam )

snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/${NN}_${TM}.snpmat.gz $NORMAL $TUMOR


snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_HeLa.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/HeLa.dd.bam &
snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_HepG2.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/HepG2.dd.bam &
snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_JurkatE6.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/JurkatE6.dd.bam &
snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_MCF7.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/MCF7.dd.bam &
snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_RAJI.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/RAJI.dd.bam &
snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_SHSY5Y.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/SHSY5Y.dd.bam &
snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_SW480.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/SW480.dd.bam &
snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_THP1.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/THP1.dd.bam &
snp-pileup -g -q15 -Q20 -P100 -r25,0 $DBSNP $OUTDIR/HG00281_HG00187.snpmat.gz bwa_map/HG00281.dd.bam bwa_map/HG00187.dd.bam &

wait ${!}



## platypus callVariants -o $OUTDIR/AllVariants.vcf --nCPU=4 --refFile=$BWA_h37 --bamFiles=$( echo $( ls bwa_map/*.dd.bam ) | tr ' ' ',')
## bgzip vcf/AllVariants.vcf && tabix vcf/AllVariants.vcf
## bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%NR:%NV]\n' vcf/AllVariants.vcf.gz > vcf/SupportingReads.txt
## cut -f 1-3,4 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/HG00281=//' | tr ":" "\t" > vcf/HG00187_HG00281.snpmat
## cut -f 1-3,5 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/HeLa.dd=//' | tr ":" "\t" > vcf/HG00187_HeLa.snpmat
## cut -f 1-3,6 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/HepG2.dd=//' | tr ":" "\t" > vcf/HG00187_HepG2.snpmat
## cut -f 1-3,7 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/JurkatE6.dd=//' | tr ":" "\t" > vcf/HG00187_JurkatE6.snpmat
## cut -f 1-3,8 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/MCF7.dd=//' | tr ":" "\t" > vcf/HG00187_MCF7.snpmat
## cut -f 1-3,9 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/RAJI.dd=//' | tr ":" "\t" > vcf/HG00187_RAJI.snpmat
## cut -f 1-3,10 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/SHSY5Y.dd=//' | tr ":" "\t" > vcf/HG00187_SHSY5Y.snpmat
## cut -f 1-3,11 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/SW480.dd=//' | tr ":" "\t" > vcf/HG00187_SW480.snpmat
## cut -f 1-3,12 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/THP1.dd=//' | tr ":" "\t" > vcf/HG00187_THP1.snpmat
## cut -f 1-2,4,3 vcf/SupportingReads.txt | sed -e 's/HG00187=//' -e 's/HG00281=//' | tr ":" "\t" > vcf/HG00281_HG00187.snpmat
## cut -f 1-2,4,5 vcf/SupportingReads.txt | sed -e 's/HG00281=//' -e 's/HeLa.dd=//' | tr ":" "\t" > vcf/HG00281_HeLa.snpmat
## cut -f 1-2,4,6 vcf/SupportingReads.txt | sed -e 's/HG00281=//' -e 's/HepG2.dd=//' | tr ":" "\t" > vcf/HG00281_HepG2.snpmat
## cut -f 1-2,4,7 vcf/SupportingReads.txt | sed -e 's/HG00281=//' -e 's/JurkatE6.dd=//' | tr ":" "\t" > vcf/HG00281_JurkatE6.snpmat
## cut -f 1-2,4,8 vcf/SupportingReads.txt | sed -e 's/HG00281=//' -e 's/MCF7.dd=//' | tr ":" "\t" > vcf/HG00281_MCF7.snpmat
## cut -f 1-2,4,9 vcf/SupportingReads.txt | sed -e 's/HG00281=//' -e 's/RAJI.dd=//' | tr ":" "\t" > vcf/HG00281_RAJI.snpmat
## cut -f 1-2,4,10 vcf/SupportingReads.txt | sed -e 's/HG00281=//' -e 's/SHSY5Y.dd=//' | tr ":" "\t" > vcf/HG00281_SHSY5Y.snpmat
## cut -f 1-2,4,11 vcf/SupportingReads.txt | sed -e 's/HG00281=//' -e 's/SW480.dd=//' | tr ":" "\t" > vcf/HG00281_SW480.snpmat
## cut -f 1-2,4,12 vcf/SupportingReads.txt | sed -e 's/HG00281=//' -e 's/THP1.dd=//' | tr ":" "\t" > vcf/HG00281_THP1.snpmat


