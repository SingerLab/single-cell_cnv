library(GenomicRanges)
library(biomaRt)
library(dplyr)

gr.cyt <- readRDS("hg19_cytoband_GRanges.rds")

### 5K bins
## get bins with gc content
k5 <- read.delim("varbin.gc.content.5k.bowtie.k50.grch37.txt")

## set up a genomic ranges object for the bins
gr.5k <- GRanges(seqnames = k5$bin.chrom, ranges = IRanges(start = k5$bin.start, end = k5$bin.end, names = 1:nrow(k5)), bin.length = k5$bin.length, gene.count = k5$gene.count, cgi.count = k5$cgi.count, dist.telomere = k5$dist.telomere, gc.content = k5$gc.content)

## add cytogenetic bands
k5 <- cbind(k5, hg19_cytoband[findOverlaps(gr.5k, gr.cyt, select = "first"), c("arm", "band", "stain")])

write.table(k5, file = "varbin.gc.content.5k.bowtie.k50.grch37.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## get ensembl genes w/hgnc
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")

grch37.genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype"), mart = grch37) %>% filter(chromosome_name %in% c(1:22, "X", "Y")) %>% arrange(chromosome_name)

grch37.genes <- GRanges(seqnames = grch37.genes$chromosome_name, ranges = IRanges(start = grch37.genes$start_position, end = grch37.genes$end_position, names = grch37.genes$ensembl_gene_id), ensembl_gene_id = grch37.genes$ensembl_gene_id, hgnc.symbol = grch37.genes$hgnc_symbol, gene.type = grch37.genes$gene_biotype)

## find and print the overlaps
( gr.5k.grch37.genes <- findOverlaps(gr.5k, grch37.genes, ignore.strand = TRUE) )
tail( gr.5k.grch37.genes <- as.data.frame(gr.5k.grch37.genes)[!duplicated(gr.5k.grch37.genes@to),]  )

any(duplicated(gr.5k.grch37.genes$subjectHits))

## assign binID to each gene
grch37.genes.5k <- grch37.genes[gr.5k.grch37.genes$subjectHits]
grch37.genes.5k$bin.id <- gr.5k.grch37.genes$queryHits
grch37.genes.5k$arm <- k5$arm[gr.5k.grch37.genes$queryHits]
grch37.genes.5k$band <- k5$band[gr.5k.grch37.genes$queryHits]
grch37.genes.5k$stain <- k5$stain[gr.5k.grch37.genes$queryHits]
grch37.genes.5k
## write it out

write.table(grch37.genes.5k, paste0("grch37.5k.gene.index.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
## sed -E 's/^chr//' grch37.5k.gene.index.txt | sort -k1V -k2n | awk '{OFS="\t"; print $9, $1,$2,$3,$6,$7,$8 }' > grch37.5k.gene.index.txt

### 20K bins
## get bins with gc content
k20 <- read.delim("varbin.gc.content.20k.bowtie.k50.grch37.txt")
## k20$bin.chrom <- paste0("chr", k20$bin.chrom)

## set up a genomic ranges object for the bins
gr.20k <- GRanges(seqnames = k20$bin.chrom, ranges = IRanges(start = k20$bin.start, end = k20$bin.end, names = 1:nrow(k20)), bin.length = k20$bin.length, gene.count = k20$gene.count, cgi.count = k20$cgi.count, dist.telomere = k20$dist.telomere, gc.content = k20$gc.content)

k20 <- cbind(k20, hg19_cytoband[findOverlaps(gr.20k, gr.cyt, select = "first"), c("arm", "band", "stain")])

write.table(k20, file = "varbin.gc.content.20k.bowtie.k50.grch37.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## find and print the overlaps
( gr.20k.grch37.genes <- findOverlaps(gr.20k, grch37.genes, ignore.strand = TRUE ) )
tail( gr.20k.grch37.genes <- as.data.frame(gr.20k.grch37.genes)[!duplicated(gr.20k.grch37.genes@to),]  )

any(duplicated(gr.20k.grch37.genes$subjectHits))

## assign binID to each gene
grch37.genes.20k <- grch37.genes[gr.20k.grch37.genes$subjectHits]
grch37.genes.20k$bin.id <- gr.20k.grch37.genes$queryHits
grch37.genes.20k$arm <- k20$arm[gr.20k.grch37.genes$queryHits]
grch37.genes.20k$band <- k20$band[gr.20k.grch37.genes$queryHits]
grch37.genes.20k$stain <- k20$stain[gr.20k.grch37.genes$queryHits]
grch37.genes.20k

## write it out
write.table(grch37.genes.20k, file = "grch37.20k.gene.index.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
## sed -E 's/^chr//' grch37.20k.gene.index.txt | sort -k1V -k2n | awk '{OFS="\t"; print $9, $1,$2,$3,$6,$7,$8 }' > grch37.20k.gene.index.txt


## 50K bins
## get bins with gc content
k50 <- read.delim("varbin.gc.content.50k.bowtie.k50.grch37.txt")

## set up a genomic ranges object for the bins
( gr.50k <- GRanges(seqnames = k50$bin.chrom, ranges = IRanges(start = k50$bin.start, end = k50$bin.end, names = 1:nrow(k50)), bin.length = k50$bin.length, gene.count = k50$gene.count, cgi.count = k50$cgi.count, dist.telomere = k50$dist.telomere, gc.content = k50$gc.content) )

## append arm band and stain
k50 <- cbind(k50, hg19_cytoband[findOverlaps(gr.50k, gr.cyt, select = "first"), c("arm", "band", "stain")])

write.table(k50, file = "varbin.gc.content.50k.bowtie.k50.grch37.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## find and print the overlaps
( gr.50k.grch37.genes <- findOverlaps(gr.50k, grch37.genes, ignore.strand = TRUE) )
tail( gr.50k.grch37.genes <- as.data.frame(gr.50k.grch37.genes)[!duplicated(gr.50k.grch37.genes@to),]  )

any(duplicated(gr.50k.grch37.genes$subjectHits))

## assign binID to each gene
grch37.genes.50k <- grch37.genes[gr.50k.grch37.genes$subjectHits]
grch37.genes.50k$bin.id <- gr.50k.grch37.genes$queryHits
grch37.genes.50k$arm <- k50$arm[gr.50k.grch37.genes$queryHits]
grch37.genes.50k$band <- k50$band[gr.50k.grch37.genes$queryHits]
grch37.genes.50k$stain <- k50$stain[gr.50k.grch37.genes$queryHits]
grch37.genes.50k


## write it out
write.table(grch37.genes.50k, file = "grch37.50k.gene.index.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
## sed -E 's/^chr//' grch37.50k.gene.index.txt | sort -k1V -k2n | awk '{OFS="\t"; print $9, $1,$2,$3,$6,$7,$8 }' > grch37.50k.gene.index.txt

