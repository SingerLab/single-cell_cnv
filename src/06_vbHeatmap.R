## plot heatmaps of varbin matrix outputs
## R.home() == "/opt/common/CentOS_7/R/R-4.0.0/lib64/R" || stop("Wrong environment, run `module load R/R-4.0.0")

args <- commandArgs(trailingOnly = TRUE)
## args <- c("--sample.name=DD4388_MR", "--io.dir=vbData/", "--fig.dir=figures/", "--aligner=bowtie")

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
      Rscript 06_vbHeatmap.R --sample.name=WD5816 --io.dir=vbData/ --fig.dir=figures/ --aligner=bowtie

")
 

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
## if(! argsL$bin.size %in% c("5k", "20k", "50k") {
##     stop("number of bins", argsL$bin.size, "does not exist")
## }

## Arg2 default
if(! dir.exists(argsL$io.dir)) {
    stop("input directory ", argsL$io.dir, " not found")
    q(save = "no")
}

## Arg3 default
if(! dir.exists(argsL$fig.dir)) {
    dir.create(argsL$fig.dir)
}


## libraries
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggplot2)
library(reshape2)
source("src/myLib.R")

## run parameters
sample.name <- argsL$sample.name
inDir  <- file.path(argsL$io.dir)
outDir <- file.path(argsL$io.dir)  ## -- will be used as input in 06_vbHeatmap.R
figDir <- file.path(argsL$fig.dir)
## bin.size <- argsL$bin.size
aligner <- argsL$aligner
bulk.pattern <- "bulk"


## read in mapd
mapd.5k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", "5k", ".k50.varbin.mapd.qc.txt", sep = "")), header = TRUE, as.is = TRUE)
mapd.5k$bin.size.kb <- 550
mapd.5k$nbins <- "5k"
mapd.20k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", "20k", ".k50.varbin.mapd.qc.txt", sep = "")), header = TRUE, as.is = TRUE)
mapd.20k$bin.size.kb <- 137
mapd.20k$nbins <- "20k"
mapd.50k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", "50k", ".k50.varbin.mapd.qc.txt", sep = "")), header = TRUE, as.is = TRUE)
mapd.50k$bin.size.kb <- 54
mapd.50k$nbins <- "50k"

## read in varbin stats
vbStats <- read.delim(file.path(inDir, paste(sample.name, "_grch37.varbin.stats.txt", sep = "")), as.is = TRUE, row.names = 1)

nn <- sapply(list(mapd.5k, mapd.20k, mapd.50k, vbStats), nrow)
names(nn) <- c("mapd.5k", "mapd.20k", "mapd.50k", "vbStats")
nn

## mapd merge
mapd.qc <- rbind(mapd.5k, mapd.20k, mapd.50k)
mapd.qc <- merge(mapd.qc, vbStats, by.x = "cellID", by.y = 0, all = TRUE)

mapd.mm <- melt(mapd.qc, id.vars = "cellID")

pdf(file.path(figDir, paste0(sample.name, "_mapd.to.bam.metrics.pdf")), width = 148/25.4, height = 105/25.4)
ggplot(mapd.qc, aes(ReadsKept, mapd, group = bin.size.kb, colour = bin.size.kb)) + geom_line() + theme_bw()
ggplot(mapd.qc, aes(DupsRemoved, mapd, group = bin.size.kb, colour = bin.size.kb)) + geom_line() + theme_bw()
ggplot(mapd.qc, aes(MedianBinCount, mapd, group = bin.size.kb, colour = bin.size.kb)) + geom_line() + theme_bw()
ggplot(mapd.qc, aes(GenomeCoverage, mapd, group = bin.size.kb, colour = bin.size.kb)) + geom_line() + theme_bw()
dev.off()


pdf(file.path(figDir, paste0(sample.name, "_mapd.stats.boxplots.pdf")), width = 148/25.4, height = 105/25.4)
## convert to ggplot boxplot + dotplot
boxplot(mapd ~ bin.size.kb, data = mapd.qc, xlab = "Bin Size (Kb)", ylab = "MAPD", cex = 0.8)
abline(h = 0.45, col = "#D40000")
text(x = 1:3, y = rep(0, 3), labels = paste("n =", round(tapply(mapd.qc$mapd, mapd.qc$bin.size.kb, function(x) sum(x < 0.45)/length(x)), 3)))
dev.off()

write.table(mapd.qc$cellID[mapd.qc$mapd > 0.45 & mapd.qc$nbins == "5k"], file = file.path(inDir, "excluded.cells.5k.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(mapd.qc$cellID[mapd.qc$mapd > 0.45 & mapd.qc$nbins == "20k"], file = file.path(inDir, "excluded.cells.20k.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(mapd.qc$cellID[mapd.qc$mapd > 0.45 & mapd.qc$nbins == "50k"], file = file.path(inDir, "excluded.cells.50k.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Read in copy number data
## 5k
cn5k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", "5k", ".k50.varbin.data.txt", sep = "")), header = TRUE, as.is = TRUE)
cn5k$chrompos <- cn5k$chrompos + 1
cn5k$abspos <- cumsum(cn5k$chrompos)

## 20k
cn20k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", "20k", ".k50.varbin.data.txt", sep = "")), header = TRUE, as.is = TRUE)
cn20k$chrompos <- cn20k$chrompos + 1
cn20k$abspos <- cumsum(cn20k$chrompos)

## 50k
cn50k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", "50k", ".k50.varbin.data.txt", sep = "")), header = TRUE, as.is = TRUE)
cn50k$chrompos <- cn50k$chrompos + 1
cn50k$abspos <- cumsum(cn50k$chrompos)

## load excluded cells
indiv <- names(cn50k[, 4:ncol(cn50k)])

## shared heatmap annotations -- chromosome color
chrLabels <- c(1:22, "X", "Y")
chrColors <- rep(c("skyblue2", "gray"), 12)
names(chrColors) <- c(1:22, "X", "Y")
chrAnno5k <- rowAnnotation(chrom = rep(c(1:22, "X", "Y"), as.numeric(table(cn5k$chrom))), col = list(chrom = chrColors), show_legend = FALSE)
chrAnno20k <- rowAnnotation(chrom = rep(c(1:22, "X", "Y"), as.numeric(table(cn20k$chrom))), col = list(chrom = chrColors), show_legend = FALSE)
chrAnno50k <- rowAnnotation(chrom = rep(c(1:22, "X", "Y"), as.numeric(table(cn50k$chrom))), col = list(chrom = chrColors), show_legend = FALSE)

## MAPD exclusion
mapd.5k$mapd.ok <- ifelse(mapd.5k$mapd <= 0.45, yes = "PASS", no = "No-pass")
mapd.20k$mapd.ok <- ifelse(mapd.20k$mapd <= 0.45, yes = "PASS", no = "No-pass")
mapd.50k$mapd.ok <- ifelse(mapd.50k$mapd <= 0.45, yes = "PASS", no = "No-pass")


## MAPD anno
hAnno5k <- HeatmapAnnotation(mapd = mapd.5k$mapd, ReadsKept = vbStats$ReadsKept/1000000, MedianBinCount  = vbStats$MedianBinCount, mapd.qc = mapd.5k$mapd.ok)
hAnno20k <- HeatmapAnnotation(mapd = mapd.20k$mapd, ReadsKept = vbStats$ReadsKept/1000000, MedianBinCount  = vbStats$MedianBinCount, mapd.qc = mapd.20k$mapd.ok)
hAnno50k <- HeatmapAnnotation(mapd = mapd.50k$mapd, ReadsKept = vbStats$ReadsKept/1000000, MedianBinCount  = vbStats$MedianBinCount, mapd.qc = mapd.50k$mapd.ok)

## shared color functions for all heatmaps
segCol <- colorRamp2(breaks = 0:5, colors = c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))
arcCol <- colorRamp2(breaks = parctan(c(0:5, 10, 20, 40, 80)), colors = c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000", "#7F7F7F", "#636363", "#515151", "#303030"))

la <- Legend(at = parctan(c(0:5, 10, 20, 40, 80)), labels = c(0:5, 10, 20, 40, ">80"), legend_gp = gpar(fill =  c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000", "#7F7F7F", "#636363", "#515151", "#303030")), title = "Integer CN")

h5k <- Heatmap(parctan(as.matrix(cn5k[, indiv])),
               cluster_rows = FALSE,
               show_column_dend = TRUE,
               column_title = paste(sample.name, "5k Heatmap"),
               show_column_names = FALSE,
               top_annotation = hAnno5k,
               left_annotation = chrAnno5k,
               row_labels = chrLabels,
               row_split = cn5k$chrom,
               row_gap = unit(0, "mm"),               
               col = arcCol,
               show_heatmap_legend = FALSE)

h20k <- Heatmap(parctan(as.matrix(cn20k[, indiv])),
               cluster_rows = FALSE,
               show_column_dend = TRUE,
               column_title = paste(sample.name, "20k Heatmap"),
               show_column_names = FALSE,
               top_annotation = hAnno20k,
               left_annotation = chrAnno20k,
               row_labels = chrLabels,
               row_split = cn20k$chrom,
               row_gap = unit(0, "mm"),               
               col = arcCol,
               show_heatmap_legend = FALSE)

h50k <- Heatmap(parctan(as.matrix(cn50k[, indiv])),
               cluster_rows = FALSE,
               show_column_dend = TRUE,
               column_title = paste(sample.name, "50k Heatmap"),
               show_column_names = FALSE,
               top_annotation = hAnno50k,
               left_annotation = chrAnno50k,
               row_labels = chrLabels,
               row_split = cn50k$chrom,
               row_gap = unit(0, "mm"),               
               col = arcCol,
               show_heatmap_legend = FALSE)

pdf(file.path(figDir, paste0(sample.name, "_vbHeatmap.quantal_%03d.pdf")), width = 298/25.4, height = 210/25.4)
draw(h5k, annotation_legend_list = packLegend(la))
draw(h20k, annotation_legend_list = packLegend(la))
draw(h50k, annotation_legend_list = packLegend(la))
dev.off()



## ggplot(mapd.mm, aes(variable, value)) + geom_boxplot()  + facet_grid(. ~ variable, scales = "free")
### attempted using umap
## p5k <- prcomp(cn5k[,4:ncol(cn5k)])
## u5k <- umap(p5k$rotation, method = "umap-learn")
## bkd5k <- dist(cn5k)
## bkd20k <- dist(cn20k)
## bkd50k <- dist(cn50k)
