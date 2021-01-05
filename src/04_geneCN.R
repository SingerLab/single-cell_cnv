## generate a gene CN for all cells
R.home() == "/opt/common/CentOS_7/R/R-4.0.0/lib64/R" || stop("Wrong environment, run `module load R/R-4.0.0`")


args <- commandArgs(trailingOnly = TRUE)
## args <- c("--sample.name=WD5816", "--io.dir=vbData/", "--fig.dir=figures/", "aligner=bowtie")

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R/QDNASeq Script
 
      Arguments:
      --sample.name=WD5816                    - character, name of sample
      --io.dir=vbData/                        - character, name of i/o directory with CN matriz
      --fig.dir=figures/                      - character, name of figurse director
      --aligner=bowtie                        - character, aligner used
      --help                                  - print this text
 
      Example:
      Rscript 04_geneCN.R --sample.name=WD5816 --io.dir=vbData/ --fig.dir=figures/ --aligner=bowtie")
 
  q(save="no")
}
##      --bin.size=(5k|20k|50k)                 - character, number of bins; for focal gene CN use 20k or 50k

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
argsL

## Arg1 default
if(! dir.exists(argsL$io.dir)) {
    stop("input directory ", argsL$io.dir, " not found")
    q(save = "no")
}

## Arg1 default
if(! dir.exists(argsL$fig.dir)) {
    dir.create(argsL$fig.dir)
}

source("src/myLib.R")

## Run parameters
sample.name <- argsL$sample.name
inDir  <- file.path(argsL$io.dir)
outDir <- file.path(argsL$io.dir)
figDir <- file.path(argsL$fig.dir)
aligner <- argsL$aligner
bulk.pattern <- "bulk"



## load copy number data
cn50k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.50k.k50.varbin.data.txt", sep = "")), header = TRUE, as.is = TRUE)

gene.index.50k <- read.delim("res/varbin/grch37.50k.gene.index.txt", as.is = TRUE)

## expanding from bin data to gene data
gene50k.hgnc <- expand2genes(cn50k[, -c(1:3)], gene.index.50k)
gene50k.ensembl <- expand2genes(cn50k[, -c(1:3)], gene.index.50k, gene.id = "ensembl_gene_id")

## save matrix in RDS form
saveRDS(gene50k.hgnc, file = file.path(outDir, paste(sample.name, "_grch37.50k.k50.varbin.hgnc.geneCN.rds", sep ="")))
saveRDS(gene50k.ensembl, file = file.path(outDir, paste(sample.name, "_grch37.50k.k50.varbin.ensembl.geneCN.rds", sep ="")))

## keep gene50k intact for future use e.g. transformation of the geneCN object
## write.table(data.frame(cellID = rownames(gene50k), gene50k), file.path(outDir, paste(sample.name, "_grch37.50k.k50.varbin.geneCN.txt", sep ="")), row.names = FALSE, quote = FALSE, sep = "\t")


## load 20k copy number data
## cn20k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.20k.k50.varbin.data.txt", sep = "")), header = TRUE, as.is = TRUE)
## gene.index.20k <- read.delim("res/varbin/grch37.20k.gencode.index.txt", as.is = TRUE)[, c(1,5)]

## gene20k <- expand2genes(cn20k, gene.index.20k)
    
## write.table(geneCN, file.path(outDir, paste(sample.name, "_grch37.20k.k50.varbin.geneCN.txt", sep ="")), row.names = FALSE, quote = TRUE, sep = "\t")
    
## load copy number data
## cn5k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.5k.k50.varbin.data.txt", sep = "")), header = TRUE, as.is = TRUE)

## 5k gene index
## gene.index.5k <- read.delim("res/varbin/grch37.5k.gene.index.txt", as.is = TRUE)

## expanding from bin data to gene data
## gene5k.hgnc <- expand2genes(cn5k[, -c(1:3)], gene.index.5k)
## gene5k.ensembl <- expand2genes(cn5k[, -c(1:3)], gene.index.5k, gene.id = "ensembl_gene_id")

## save matrix in RDS form
## saveRDS(gene5k.hgnc, file = file.path(outDir, paste(sample.name, "_grch37.5k.k50.varbin.hgnc.geneCN.rds", sep ="")))
## saveRDS(gene5k.ensembl, file = file.path(outDir, paste(sample.name, "_grch37.5k.k50.varbin.ensembl.geneCN.rds", sep ="")))


## geneCNT <- parctan(geneCN)
## geneDist <- dist(t(geneCNT))
## geneClust <- hclust(geneDist)
## sampleDist <- dist(geneCNT)
## sampleClust <- hclust(sampleDist)
## Heatmap(geneCNT)
