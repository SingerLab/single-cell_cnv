## plot heatmaps of CN5K
library(ggplot2)
library(gplots)
library(vegan)
library(copynumber)
library(pheatmap)
library(multcomp)
library(viridis)
library(colourvalues)
library(ape)
library(dendsort)
source("src/myLib.R")

## Run parameters
sample.name <- "WD9048R_MR_T"
inDir  <- file.path("vbData")
outDir <- file.path("figures")
mcpDir <- file.path("multComp")
figDir <- file.path("figures")
bin.size <- "5k"
aligner <- "bowtie"

dir.exists(figDir) || dir.create(figDir)
dir.exists(mcpDir) || dir.create(mcpDir)
dir.exists(outDir) || dir.create(outDir)

cn5k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.nobad.varbin.lowratio.data.txt", sep = "")), header = TRUE, as.is = TRUE)
cn5k$chrompos <- cn5k$chrompos + 1
cn5k$abspos <- cumsum(cn5k$chrompos)

cn20k <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", "20k", ".k50.nobad.varbin.lowratio.data.txt", sep = "")), header = TRUE, as.is = TRUE)
cn20k$chrompos <- cn20k$chrompos + 1
cn20k$abspos <- cumsum(cn20k$chrompos)

## load excluded dna
## exclBulk <- c(readLines(file.path(inDir, "excluded.bulk.txt")))
## cn5k <- cn5k[, !names(cn5k) %in% exclBulk]

## Bulk names
Bulk <-  names(cn5k[, 4:ncol(cn5k)])

## bulk annotations -- best if done manually
## include rownames in the first colum and minimum required columns are :
##(rownames) 	bulkID	subsample	plate	subtype	barcode	path.comments
bulkAnnot <- read.delim(file.path(inDir, paste(sample.name, "_bulkAnnotations.txt", sep = "")), header = TRUE, sep = "\t", as.is = TRUE, row.names = 1)

subtypeColor <- RColorBrewer::brewer.pal(length(unique(bulkAnnot$subtype)), "Set3")
names(subtypeColor) <- unique(bulkAnnot$subtype)
subsampleColor <- rainbow(nrow(bulkAnnot))
names(subsampleColor) <- bulkAnnot$bulkID
plateColor <- sample(colors(distinct = FALSE), length(unique(bulkAnnot$plate)))
names(plateColor) <- unique(bulkAnnot$plate)

## importing geneCN
bulkCN <- read.delim(file.path(inDir, paste(sample.name, "_grch37.50k.k50.nobad.varbin.gene_focal_CN.txt", sep = "")), stringsAsFactors = FALSE)
## annotating bulk type
bulkCN$bulk.type <- "normal.bulk"
bulkCN$bulk.type[bulkCN$CDK4 > 2.5 | bulkCN$CDK4 < 1.5] <- "tumor.bulk"
bulkCN$bulk.type[bulkCN$MDM2 > 2.5 | bulkCN$MDM2 < 1.5] <- "tumor.bulk"

## importing FGA data
bulkFGA <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.nobad.varbin.FGA.txt", sep = "")), stringsAsFactors = FALSE)

## importing ploidy
bulkPloidy <- read.delim(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.nobad.varbin.quantal.stats.txt", sep = "")), stringsAsFactors = FALSE)
( names(bulkPloidy)[2:4] <- paste(names(bulkPloidy)[2:4], "5k", sep = ".") )

## merging to geneCN
bulkCN <- merge(bulkAnnot, bulkCN, by.x = "bulkID", by.y = "cellID")
bulkCN <- merge(bulkCN, bulkFGA, by.x = "bulkID", by.y = "cellID")
bulkCN <- merge(bulkCN, bulkPloidy, by.x = "bulkID", by.y = "cellID")
rownames(bulkCN) <- bulkCN$bulkID

## assuming NA are 0
if(any(is.na(bulkCN))) {
    bulkCN$n.segments.altered[is.na(bulkCN$n.segments.altered)] <- 0
    bulkCN$bp.altered[is.na(bulkCN$bp.altered)] <- 0
    bulkCN$n.segments.altered.20mb[is.na(bulkCN$n.segments.altered.20mb)] <- 0
    bulkCN$fga[is.na(bulkCN$fga)] <- 0
}

## removing rows with missing values
## bulkCN <- bulkCN[rowSums(is.na(bulkCN)) == 0, ]

## asigning colors to t-cells
bulkCN$bulkTypeColor <- colour_values(bulkCN$bulk.type, palette = "cividis")

bulkTypeColor <- unique(bulkCN$bulkTypeColor)
names(bulkTypeColor) <- unique(bulkCN$bulk.type)

## asigning colors gene and fga  -- limited to 41 and 0.5
tmp1 <- bulkCN$MDM2
tmp1[tmp1 >= 41] <- 41
bulkCN$MDM2.color <- colour_values(tmp1, palette = "magma")

tmp2 <- bulkCN$CDK4
tmp2[tmp2 >= 40] <- 41
bulkCN$CDK4.color <- colour_values(tmp2, palette = "magma")

tmp3 <- bulkCN$fga
tmp3[tmp3 >= 0.35] <- 0.35
bulkCN$fgaColor <- colour_values(tmp3, palette = "viridis")

rm(tmp1, tmp2, tmp3)

## genome mapping annotation
chrBreaks <- cumsum(table(cn5k$chrom))
chrCol <- rep(c("#9E9E9E", "#64B5F6"), 12)
names(chrCol) <- c(1:22, "X", "Y")
chrCol <- chrCol[cn5k$chrom]
chrLabels <- rep(NA, nrow(cn5k))
chrLabels[c(1, chrBreaks[1:23])+20] <- c(1:22, "X", "Y")

cnr <- as.matrix(cn5k[, c(Bulk)])

cnm <- log2(cnr)

## 20k
cnr20k <- as.matrix(cn20k[, c(Bulk)])

## heatmap## heatmap
## color choices --
## #F2D200 - yellow 62 - homozygous deletion
## #318CE7 - french racing blue - cn loss
## #FFFFFF - white - 2n -- check perl white #EAE0C8	
## #FF523F - ducati corse orange - gain
## #D40000 - rosso corse - amplification
## #032ADD - another type of blue or yellow, can't remember :)
## "#7F7F7F", "#636363", "#515151", "#303030" - gray scale

## calculating distance of all bulks
cdd <- dist(t(cnm))
mdd <- cmdscale(cdd)
hcdd <- hclust(cdd)

## deciding where to cut
## number of clusters by tree height
tree.height <- 1:24
lapply(tree.height, function(he) {
    cct <- cutree(hcdd, h = he)
##    summaryp(cct)
    as.data.frame(table(cct))
    })

## cut tree screening 
png(file.path(figDir, paste(sample.name, "_bulk_%03d.", aligner, bin.size, "_clusterOptimization.png", sep = "")), width = 410, height = 298, units = "mm", res = 400)
tree.height = 1:16
sapply(tree.height, function(he) {
    cct <- cutree(hcdd, h = he)
    clColor <- rainbow(length(unique(cct)))
    names(clColor) <- paste("B_C", unique(cct), sep = "")
    myClusterSideBar <- clColor[cct]
    heatmap.2(cnm[, c(Bulk)], Rowv = FALSE, Colv = hcdd,
              dendrogram = "col", trace = "none",
              col = colorRampPalette(viridis(n = 10)),
              main = paste(sample.name, " -- WDLS (n =", length(Bulk), "); tree height =", he),  yaxt = "n",
              ColSideColors = myClusterSideBar, cexCol = 0.8, 
              offsetCol = 1,
              RowSideColors = chrCol, cexRow = 1.4,
              labRow = chrLabels, offsetRow = (104 * -1),
              keysize = 0.6)
    legend("right", legend = names(clColor), col = clColor, pch = 15, cex = .6, ncol = 2)
})
dev.off()


## cutting the tree
tree.height <- 11
cn5k.clust <- cutree(hcdd, h = tree.height)
bulkCN$cluster <- paste("B_C", cn5k.clust[rownames(bulkCN)], sep = "")
bulkCN$cluster[bulkCN$cluster == "B_CNA"] <- NA
bulkCN$cluster <- factor(bulkCN$cluster)

## use civids instead of rainbow
clColor <- cividis(length(unique(bulkCN$cluster)))
names(clColor) <- sort(unique(bulkCN$cluster))
bulkCN$clusterColor <- clColor[bulkCN$cluster]
bulkCN$subtypeColor <- subtypeColor[bulkCN$subtype]
bulkCN$subsampleColor <- subsampleColor[gsub("([0-9]).2$", "\\1", bulkCN$subsample)]

nodePar  <- list(list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = bulkCN$clusterColor[hcdd$lables]))
dhcdd <- as.dendrogram(hcdd)
                
pdf(file.path(figDir, paste(sample.name, "_MR_bulk_dendrogram.pdf", sep = "")), width = 298/25.4, height = 149/25.4)
plot(hcdd, cex = 0.1, xlab = "", sub = "", lwd = 0.8)
plot(as.phylo(hcdd), type = "unrooted", cex = 0.4, tip.color = bulkCN[hcdd$labels, "clusterColor"], lwd = 1.3)
dev.off()


## check multiple comparisons directory
dir.exists(mcpDir) || dir.create(mcpDir)

## glm all comparisons
glmMDM2 <- glm(MDM2 ~ cluster, data = bulkCN)
MDM2.mc <- glht(glmMDM2, linfct = mcp(cluster = "Tukey"))
MDM2.smry <- summary(MDM2.mc)
write.table(table.glht(MDM2.mc), file = file.path(mcpDir, paste0(sample.name, "_MDM2_multiple_comparisons.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## glm all comparisons
glmCDK4 <- glm(CDK4 ~ cluster, data = bulkCN)
CDK4.mc <- glht(glmCDK4, linfct = mcp(cluster = "Tukey"))
CDK4.smry <- summary(CDK4.mc)
write.table(table.glht(CDK4.mc), file = file.path(mcpDir, paste0(sample.name, "_CDK4_multiple_comparisons.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## glm all comparisons
glmFGA <- glm(fga ~ cluster, data = bulkCN)
FGA.mc <- glht(glmFGA, linfct = mcp(cluster = "Tukey"))
FGA.smry <- summary(FGA.mc)
write.table(table.glht(FGA.mc), file = file.path(mcpDir, paste0(sample.name,"_FGA_multiple_comparisons.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


pdf(file.path(figDir, paste0(sample.name, "_Dosage_Cluster_boxplot.pdf")), width = 149/25.4, height = 105/25.4)
ggplot(bulkCN, aes(cluster, MDM2)) + geom_boxplot(outlier.color = "blue", outlier.size = 0.5) +
    geom_dotplot(alpha = 0.6, binaxis = "y", stackdir="center", dotsize = 0.2) +
    stat_summary(fun.data = stat_box_data, geom = "text") +
    xlab("Cluster ID") + ylab("MDM2 CN") +
    labs(title = "MDM2 copy number state") +
    theme_bw()
ggplot(bulkCN, aes(cluster, CDK4)) + geom_boxplot(outlier.color = "blue", outlier.size = 0.5) +
    geom_dotplot(alpha = 0.6, binaxis = "y", stackdir="center", dotsize = 0.2) +
    stat_summary(fun.data = stat_box_data, geom = "text") +
    xlab("Cluster ID") + ylab("CDK4 CN") +
    labs(title = "CDK4 copy number state") +
    theme_bw()
ggplot(bulkCN, aes(cluster, fga)) + geom_boxplot(outlier.color = "blue", outlier.size = 0.5) +
    geom_dotplot(alpha = 0.6, binaxis = "y", stackdir="center", dotsize = 0.2, drop = TRUE) +
    xlab("Cluster ID") + ylab("FGA") +
    labs(title = "Fraction of the Genome Altered (CN != 2n)") +
    theme_bw()
ggplot(bulkCN, aes(cluster, fga)) + geom_boxplot(outlier.color = "blue", outlier.size = 0.5) +
    geom_dotplot(alpha = 0.6, binaxis = "y", stackdir="center", dotsize = 0.2, drop = TRUE) +
    xlab("Cluster ID") + ylab("FGA") + ylim(c(0,0.5)) +
    labs(title = "Fraction of the Genome Altered (CN != 2n)") +
    theme_bw()
ggplot(bulkCN, aes(cluster, ploidy.5k)) +
    geom_boxplot(outlier.color = "blue", outlier.size = 0.5) +
    geom_jitter(alpha = 3/10, width = .2, size = 0.4) +
    xlab("Cluster ID") + ylab("Estiamted Ploidy") +
    labs(title = "Estimated bulk ploidy") +
    theme_bw()
dev.off()


pdf(file.path(figDir, paste0(sample.name,  "_Dosage_Cluster_cor.pdf")))
pairs(bulkCN[bulkCN$bulk.type == "tumor.bulk", c("fga", "MDM2", "CDK4")],
      panel = panel.lm, col = adjustcolor(4, .3), pch = 19,
      cex = 0.8, font.labels = 1.8, lower.panel = panel.cor)
plot(MDM2 ~ CDK4, data = bulkCN, col = bulkCN$bulkTypeColor, xlim = c(0,100), ylim = c(0,100))
abline(a = 0, b = 1)
abline(a = lm(MDM2 ~ CDK4, data = bulkCN[bulkCN$bulk.type == "tumor.bulk", ]), col = "#D40000")
legend("topleft", legend = names(bulkTypeColor), col = bulkTypeColor, pch = 1, pt.cex = 2)
dev.off()


## plotting heatmaps
png(file.path(figDir, paste(sample.name, "_MR_bulk.", aligner, bin.size, "_dosageClusters_%03d.png", sep = "")), width = 410, height = 298, units = "mm", res = 400)
heatmap.2(cnm[, bulkCN$bulkID], Rowv = FALSE, Colv = TRUE,
          dendrogram = "col", trace = "none",
          col = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
          main = paste(sample.name, " -- WDLS Bulk DNA (n =", sum(bulkCN$bulkID %in% Bulk), "); tree height =", tree.height),  yaxt = "n",
          ColSideColors = bulkCN$clusterColor, cexCol = 0.8, 
          labCol = bulkCN$subsample, offsetCol = 1,
          RowSideColors = chrCol, cexRow = 1.4,
          labRow = chrLabels, offsetRow = (122 * -1),
          keysize = 0.6)
legend("right", legend = names(clColor), col = clColor, pch = 15, cex = 0.8, ncol = 2)
heatmap.2(cnm[, bulkCN$bulkID], Rowv = FALSE, Colv = TRUE,
          dendrogram = "col", trace = "none",
          col = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
          main = paste(sample.name, " -- WDLS Bulk DNA (n =", sum(bulkCN$bulkID %in% Bulk), "); tree height =", tree.height),  yaxt = "n",
          ColSideColors = bulkCN$clusterColor, cexCol = 1.1, 
          labCol = gsub("WD4495[TN_\\.]", "", bulkCN$subsample), offsetCol = -102,
          srtCol=0,   adjCol = c(0.5,1),
          RowSideColors = chrCol, cexRow = 1.4,
          labRow = chrLabels, offsetRow = (122 * -1),
          keysize = 0.6)
legend("right", legend = names(clColor), col = clColor, pch = 15, cex = 0.8, ncol = 2)
dev.off()


png(file.path(figDir, paste(sample.name, "_MR_bulk.", aligner, bin.size, "_heatmap_%03d.png", sep = "")), width = 410, height = 298, units = "mm", res = 400)
heatmap.2(cnm[, bulkCN$bulkID], Rowv = FALSE, Colv = TRUE,
          dendrogram = "col", trace = "none",
          col = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
          main = paste(sample.name, " -- WDLS Bulk DNA"),  yaxt = "n",
          ColSideColors = bulkCN$subtypeColor, cexCol = 0.8, 
          offsetCol = 1,
          RowSideColors = chrCol, cexRow = 1.4,
          labRow = chrLabels, offsetRow = (124 * -1),
          keysize = 0.6)
legend("right", legend = names(subtypeColor), col = subtypeColor, pch = 15, cex = 1, ncol = TRUE)
heatmap.2(cnm[, bulkCN$bulkID], Rowv = FALSE, Colv = TRUE,
          dendrogram = "col", trace = "none",
          col = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
          main = paste(sample.name, " -- WDLS Bulk DNA"),  yaxt = "n",
          ColSideColors = bulkCN$subsampleColor, cexCol = 0.8, 
          offsetCol = 1,
          RowSideColors = chrCol, cexRow = 1.4,
          labRow = chrLabels, offsetRow = (124 * -1),
          keysize = 0.6)
legend("right", legend = names(subsampleColor), col = subsampleColor, pch = 15, cex = 0.8, ncol = 2)
heatmap.2(cnm[, bulkCN$bulkID], Rowv = FALSE, Colv = TRUE,
          dendrogram = "col", trace = "none",
          col = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
          main = paste(sample.name, " -- WDLS Bulk DNA CDK4 CN"),  yaxt = "n",
          ColSideColors = bulkCN$CDK4.color, cexCol = 0.8, 
          offsetCol = 1,
          RowSideColors = chrCol, cexRow = 1.4,
          labRow = chrLabels, offsetRow = (124 * -1),
          keysize = 0.6)
legend("topright", legend = c(bulkCN[which.min(bulkCN$CDK4), "CDK4"], bulkCN[which.max(bulkCN$CDK4), "CDK4"]), col = c(bulkCN[which.min(bulkCN$CDK4), "CDK4.color"], bulkCN[which.max(bulkCN$CDK4), "CDK4.color"]), pch = 15, cex = .8)
heatmap.2(cnm[, bulkCN$bulkID], Rowv = FALSE, Colv = TRUE,
          dendrogram = "col", trace = "none",
          col = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
          main = paste(sample.name, " -- WDLS Bulk DNA MDM2 CN"),  yaxt = "n",
          ColSideColors = bulkCN$MDM2.color, cexCol = 0.8, 
          offsetCol = 1,
          RowSideColors = chrCol, cexRow = 1.4,
          labRow = chrLabels, offsetRow = (124 * -1),
          keysize = 0.6)
legend("topright", legend = c(bulkCN[which.min(bulkCN$MDM2), "MDM2"], bulkCN[which.max(bulkCN$MDM2), "MDM2"]), col = c(bulkCN[which.min(bulkCN$MDM2), "MDM2.color"], bulkCN[which.max(bulkCN$MDM2), "MDM2.color"]), pch = 15, cex = .8)
heatmap.2(cnm[, bulkCN$bulkID], Rowv = FALSE, Colv = TRUE,
          dendrogram = "col", trace = "none",
          col = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
          main = paste(sample.name, " -- WDLS Bulk DNA FGA"),  yaxt = "n",
          ColSideColors = bulkCN$fgaColor, cexCol = 0.8, 
          offsetCol = 1,
          RowSideColors = chrCol, cexRow = 1.4,
          labRow = chrLabels, offsetRow = (124 * -1),
          keysize = 0.6)
legend("topright", legend = c(bulkCN[which.min(bulkCN$fga), "fga"], bulkCN[which.max(bulkCN$fga), "fga"]), col = c(bulkCN[which.min(bulkCN$fga), "fgaColor"], bulkCN[which.max(bulkCN$fga), "fgaColor"]), pch = 15, cex = 0.8)
heatmap.2(cnm[, bulkCN$bulkID], Rowv = FALSE, Colv = TRUE,
          dendrogram = "col", trace = "none",
          col = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
          main = paste(sample.name, " -- WDLS Bulk DNA FGA"),  yaxt = "n",
          ColSideColors = bulkCN$bulkTypeColor, cexCol = 0.8, 
          offsetCol = 1,
          RowSideColors = chrCol, cexRow = 1.4,
          labRow = chrLabels, offsetRow = (124 * -1),
          keysize = 0.6)
legend("topright", legend = names(bulkTypeColor), col = bulkTypeColor, pch = 15, cex = 0.8)
heatmap.2(t(as.matrix(bulkCN[, c("MDM2", "CDK4")])),
          trace = "none",
          col = viridis(n = nrow(bulkCN), option = "D"),
          main = paste(sample.name, " -- WDLS Bulk DNA gene CN"),  yaxt = "n",
          ColSideColors = bulkCN[, "clusterColor"])
heatmap.2(t(as.matrix(bulkCN[, c("MDM2", "CDK4")])),
          trace = "none",
          col = viridis(n = nrow(bulkCN), option = "A"),
          main = paste(sample.name, " -- WDLS Bulk DNA gene CN"),  yaxt = "n",
          ColSideColors = bulkCN[, "clusterColor"])
dev.off()



## pheatmap annotation table
heatAnnot <- bulkCN[, c("cluster", "CDK4", "MDM2", "fga", "n.segments.20mb.deletion", "ploidy.5k", "bulk.type", "subtype", "subsample")]
names(heatAnnot)[9] <- "region"
heatAnnot <- droplevels(heatAnnot)

## pheatmap annotation colors -- does not work; use default colors
heatAnnotColors <- data.frame(cluster = rep(NA, times= nrow(bulkCN)))
heatAnnotColors$cluster <- colour_values(bulkCN$cluster, palette = "cividis")
heatAnnotColors$CDK4 <- colour_values(bulkCN$CDK4, palette = "magma")
heatAnnotColors$MDM2 <- colour_values(bulkCN$MDM2, palette = "magma")
heatAnnotColors$fga <- colour_values(bulkCN$fga)
heatAnnotColors$n.segments.20mb.deletion <- colour_values(bulkCN$n.segments.20mb.deletion)
heatAnnotColors$ploidy.5k <- colour_values(bulkCN$ploidy.5k, palette = "blues")
heatAnnotColors$bulk.type <- colour_values(bulkCN$bulk.type, palette = "cividis")
heatAnnotColors$subtype <- colour_values(bulkCN$subtype, palette = "diverge_hsv")
heatAnnotColors$region <- colour_values(bulkCN$subsample, palette = "diverge_hsv")
heatAnnotColors <- droplevels(heatAnnotColors)

dosageClust <- hclust(dist(t(cnm[,bulkCN$bulkID])))
sdClust <- as.hclust(dendsort(as.dendrogram(dosageClust), isReverse = TRUE))

chrAnnot <- data.frame(chrom = factor(gsub("24", "Y", gsub("23", "X", cn5k$chrom)), levels = c(1:22, "X", "Y")), row.names = cn5k$abspos)
chrAnnotCol <- data.frame(chrCol, row.names = cn5k$abspos)
rownames(cnm) <- cn5k$abspos

png(file.path(figDir, paste0(sample.name, "_MR_bulk.%03d_preetyHeatmap.png")), width = 410, height = 330, units = "mm", res = 400)
pheatmap(mat = cnm[, Bulk],
         color = colorRampPalette(c("#F2D200", "#318CE7", "#FFFFFF", "#FFA89F", "#FF523F", "#D40000"))(6),
         main = paste(sample.name, "-- WDLS Bulk DNA (n =", nrow(bulkCN), "); tree height =", tree.height),
         border_color = NA, drop_levels = TRUE,
         cluster_cols = dosageClust, cluster_rows = FALSE,
         clustering_method = "complete",
         annotation_col = heatAnnot,
         labels_row = chrLabels,
##         annotation_color = heatAnnotColors,
         annotation_row = chrAnnot,
         show_rownames = FALSE,
         fontsize_col = 2)
dev.off()


save.image("06_bulkHeatmap.rda")
