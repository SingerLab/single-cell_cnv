## assesing tcell infiltration
## libraries
library(colourvalues)

## run parameters
sample.name <- "WD5816"
inDir  <- file.path("vbData/")
outDir <- file.path("tcell/")  ## -- will be used as input in 06_vbHeatmap.R
bin.size <- "50k"
aligner <- "bowtie"

## read in matrix data -- re-bulild seg file to have copy number data instead of ratios
## cn50k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.nobad.varbin.data.txt", sep = "")), header = TRUE, as.is = TRUE)
## seg50k <- read.table(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.nobad.varbin.short.cbs.seg", sep = "")), header = FALSE, col.names = c("cellID", "chrom", "start", "end", "nmark", "ratio"), as.is = TRUE)

## load excluded cells
## exclCells <- c(readLines("excluded.cells.txt"), readLines("aberrant.cells.txt"))
## cn50k <- cn50k[, !names(cn50k) %in% exclCells]

## get cell names
## Cells <- names(cn50k[, 4:ncol(cn50k)])
## cnr <- round(cn50k)
## cnr <- cnr[, -c(3)]

tcell <- readLines(file.path(inDir, "putative.R9.50k_2x.t.cell.txt"))

## cell annotations
## source("src/06_cellAnnot.R", verbose = TRUE)
cellAnnot <- read.delim(file.path(inDir, paste(sample.name, "_cellAnnotations.txt", sep = "")), header = TRUE, sep = "\t", as.is = TRUE, row.names = 1)
## cellAnnot <- cellAnnot[cellAnnot$cellID %in% c(Cells),]

subtypeColor <- readRDS(file.path(inDir, "subtypeColor.rds"))
subsampleColor <- readRDS(file.path(inDir, "subsampleColor.rds"))
plateColor <- readRDS(file.path(inDir, "plateColor.rds"))

## importing geneCN
geneCN <- read.delim(file.path(inDir, paste(sample.name, "_grch37.50k.k50.nobad.varbin.gene_focal_CN.txt", sep = "")), stringsAsFactors = FALSE)
geneCN$cell.type <- "normal.cell"
geneCN$cell.type[geneCN$CDK4 > 2.5 | geneCN$CDK4 < 1.5] <- "tumor.cell"
geneCN$cell.type[geneCN$MDM2 > 2.5 | geneCN$MDM2 < 1.5] <- "tumor.cell"

## importing FGA data
cellFGA <- read.delim(file.path(inDir, paste(sample.name, "_grch37.5k.50k.nobad.varbin.FGA.txt", sep = "")), stringsAsFactors = FALSE)

## merging to geneCN
geneCN <- merge(cellAnnot, geneCN, by = "cellID", all = TRUE)
geneCN <- merge(geneCN, cellFGA, by = "cellID", all = TRUE)
rownames(geneCN) <- geneCN$cellID

## assuming NA are 0
geneCN$n.segments.altered[is.na(geneCN$n.segments.altered)] <- 0
geneCN$bp.altered[is.na(geneCN$bp.altered)] <- 0
geneCN$n.segments.altered.20mb[is.na(geneCN$n.segments.altered.20mb)] <- 0
geneCN$fga[is.na(geneCN$fga)] <- 0

## removing rows with missing values
## geneCN <- geneCN[rowSums(is.na(geneCN)) == 0, ]

## marking t-cells
tcell <- readLines(file.path(inDir, "putative.R9.50k_2x.t.cell.txt"))
geneCN$cell.type[geneCN$cellID %in% tcell] <- "t-cell"

## asigning colors to t-cells
geneCN$cellTypeColor <- colour_values(geneCN$cell.type, palette = "cividis")

cellTypeColor <- unique(geneCN$cellTypeColor)
names(cellTypeColor) <- unique(geneCN$cell.type)

## asigning colors gene and fga  -- limited to 41 and 0.5
tmp1 <- geneCN$MDM2
tmp1[tmp1 >= 41] <- 41
geneCN$MDM2.color <- colour_values(tmp1, palette = "magma")

tmp2 <- geneCN$CDK4
tmp2[tmp2 >= 40] <- 41
geneCN$CDK4.color <- colour_values(tmp2, palette = "magma")

tmp3 <- geneCN$fga
tmp3[tmp3 >= 0.5] <- 0.5
geneCN$fgaColor <- colour_values(tmp3, palette = "viridis")

rm(tmp1, tmp2, tmp3)


## used for the CSHL abstract -- 
## tcell.ct <- matrix(c(75,97,94,531), nrow = 2, dimnames = list("cell.type" = c("T-Cell", "Non-T-Cell"), "sample" = c("DDLS", "DDLS c/w RMS")))
tcell <- geneCN[, c("cell.type", "subtype")]

tcell.ct <- table(tcell$cell.type, tcell$subtype)[, 1:3]
( tcell.fexact <- fisher.test(tcell.ct, simulate.p.value = TRUE, B=10000) )


pdf(file.path(outDir, "tcell_barplot.pdf"))
barplot(tcell.ct, col = cellTypeColor[rownames(tcell.ct)], legend.text = TRUE, xlab = "Histological Subtype", ylab = "Number of Cells")
## legend("topleft", legend = names(cellTypeColor), col = cellTypeColor, pch = 15, pt.cex = 3)
dev.off()


## 

## summary(aov(ploidy.5k ~ cell.type, data = cellAnnot))

## pdf("T_N_ploidy.pdf")
## ggplot(cellAnnot, aes(cell.type, ploidy.5k)) + geom_boxplot(fill = "#DFDFDF", outlier.colour = "#D40000") +
##     geom_jitter(alpha = 5/10, width = 0.18) +
##     theme_minimal()
## dev.off()


## pdf("12q_Plate.integerCN.violin.pdf")
## ggplot(geneCNm[geneCNm$cell.type == "tumor.cell",], aes(plate, copy.number)) + geom_violin(alpha = 3/10) +
##     geom_dotplot(alpha = 6/10, binaxis = 'y', stackdir = 'center', dotsize = 0.08) +
##     facet_wrap(~ gene, scales = "free_y") + theme_linedraw()
## dev.off()

