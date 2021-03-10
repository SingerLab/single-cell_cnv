## estimating the fraction of genome altered and frequency per cells
## R.home() == "/opt/common/CentOS_7/R/R-4.0.0/lib64/R" || stop("Wrong environment, run `module load R/R-4.0.0`")

args <- commandArgs(trailingOnly = TRUE)
## args <- c("--sample.name=WD5816", "--bin.size=5k", "--io.idr=vbData/", "--fig.dir=figures/", "--aligner=bowtie")

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
      --bin.size=(5k|20k|50k)                 - character, number of bins, for fga use 5k or 20k
      --io.dir=vbData/                        - character, name of i/o directory with CN matrix
      --fig.dir=figures/                      - character, name of figurse director
      --aligner=bowtie                        - character, aligner used
      --help                                  - print this text
 
      Example:
      Rscript 04_geneCN.R --sample.name=WD5816 --bin.size=5k --io.dir=vbData/ --fig.dir=figures/ --aligner=bowtie")
 
  q(save="no")
}


## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
argsL

## Arg1 default
if(! argsL$bin.size %in% c("5k", "20k", "50k")) {
    stop("number of bins", argsL$bin.size, "does not exist")
}

## Arg2 default
if(! dir.exists(argsL$io.dir)) {
    stop("input directory ", argsL$io.dir, " not found")
    q(save = "no")
}

## Arg3 default
if(! dir.exists(argsL$fig.dir)) {
    dir.create(argsL$fig.dir)
}

## libraries and resources
library(copynumber)
source("src/myLib.R")

## run parameters
sample.name <- argsL$sample.name
inDir  <- file.path(argsL$io.dir)
outDir <- file.path(argsL$io.dir)  ## -- will be used as input in 06_vbHeatmap.R
figDir <- file.path(argsL$fig.dir)
bin.size <- argsL$bin.size
aligner <- argsL$aligner
bulk.pattern <- "bulk"

dir.exists(figDir) || dir.create(figDir)

## read in matrix data -- re-bulild seg file to have copy number data instead of ratios
cn5k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.varbin.data.txt", sep = "")), header = TRUE, as.is = TRUE)
lr5k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.varbin.lowratio.data.txt", sep = "")), header = TRUE, as.is = TRUE)
## seg5k <- read.table(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.varbin.short.cbs.seg", sep = "")), header = FALSE, col.names = c("cellID", "chrom", "start", "end", "nmark", "ratio"), as.is = TRUE)

## load excluded cells
## exclCells <- c(readLines("excluded.cells.txt"), readLines("aberrant.cells.txt"))
## cn5k <- cn5k[, !names(cn5k) %in% exclCells]

## get cell names
Cells <- names(cn5k[, 4:ncol(cn5k)])
cnr <- round(cn5k)
cnr <- cnr[, -c(3)]

Bulks <- names(lr5k[, grep(bulk.pattern, names(lr5k))])
## bkr <- lr5k[, c("chrom", "chrompos", Bulks)]
## bkr[,Bulks] <- log2(lr5k[, Bulks])

## convert matrix into seg data using copynumber
## wins <- winsorize(data = cnr, outliers = TRUE, verbose = FALSE)
seg5k <- pcf(data = cnr, gamma = 12, verbose = FALSE)
seg5k$mean <- round(seg5k$mean)
seg5k <- seg5k[order(seg5k$sampleID),]
seg5k$segment.length <- seg5k$end.pos - seg5k$start.pos

## defining altered segments  -- not diploid, not chrY
altered5k <- seg5k[seg5k$mean != 2 & seg5k$chrom != 24,]

## generating alteration summary based on the number of altered segments, and number of bases altered
alterationSummary <- data.frame(table(altered5k$sampleID))
names(alterationSummary) <- c("cellID", "n.segments.altered")
alterationSummary$bp.altered <- tapply(altered5k$segment.length, altered5k$sampleID, sum)

## bringing non-altered cells to have the complete set of cells in FGA
if(sum(!unique(seg5k$sampleID) %in% unique(altered5k$sampleID)) > 0) {
    notAltered5k <- data.frame(cellID = unique(seg5k$sampleID)[!unique(seg5k$sampleID) %in% unique(altered5k$sampleID)], n.segments.altered = 0, bp.altered = 0)
} else {
    notAltered5k <- data.frame(cellID = character(), n.segments.altered = integer(), bp.altered = integer())
}

## binding the altered with non-altered cells
alterationSummary <- rbind(alterationSummary, notAltered5k)

## looking into cells with alterations greater than 20 Mb
mb20 <- as.data.frame(table(seg5k[seg5k$mean != 2 & seg5k$chrom != 24 & seg5k$segment.length >= 20000000, "sampleID"]))
if(nrow(mb20) != 0)  names(mb20) <- c("cellID", "n.segments.altered.20mb")

## looking at large deletions (> 20 Mb)
mb20L <- as.data.frame(table(seg5k[seg5k$mean == 0 & seg5k$chrom != 24 & seg5k$segment.length >= 20000000, "sampleID"]))
if(nrow(mb20L) != 0) names(mb20L) <- c("cellID", "n.segments.20mb.deletion")

## merging to alteration summary
if(nrow(mb20) != 0)  alterationSummary <- merge(alterationSummary, mb20, by = "cellID", all = TRUE)
if(nrow(mb20L) != 0) alterationSummary <- merge(alterationSummary, mb20L, by = "cellID", all = TRUE)

## test this first
if(!"n.segments.20mb.deletion" %in% names(alterationSummary)) alterationSummary$n.segments.20mb.deletion <- 0

## FGA = Fraction of Genome Altered / proportion of the genome not in diploid state
## usisng goldenPath genome length
alterationSummary$fga <- alterationSummary$bp.altered / 3098825702

alterationSummary[is.na(alterationSummary)] <- 0

write.table(alterationSummary, file = file.path(outDir, paste(sample.name, "_grch37.", bin.size, ".k50.varbin.FGA.txt", sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)

## write.table(alterationSummary[grep(bulk.pattern, alterationSummary$cellID), ], file = file.path(outDir, paste(sample.name, "_bulk_grch37.", bin.size, ".k50.varbin.FGA.txt", sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)
## head(alterationSummary)

pdf(file.path(figDir, paste(sample.name, "_grch37.", bin.size, ".k50.cna.pdf", sep = "")), width = 298/25.4, height = 210/25.4)
plotAberration(seg5k, thres.gain = 0.1, thres.loss = -0.2)
dev.off()

seg5km <- multipcf(data = cnr, verbose = FALSE, gamma = 12)
pdf(file.path(figDir, paste(sample.name, "_grch37.", bin.size, ".k50.cna.pdf", sep = "")), width = 210/25.4, height = 210/25.4)
plotCircle(seg5km, thres.gain = 0.15)
dev.off()
