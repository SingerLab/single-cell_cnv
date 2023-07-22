## Table for export results of multiple comparison (post hoc Tukey)
## table.glht -- Source: Modified from https://gist.github.com/cheuerde/3acc1879dc397a1adfb0
## https://gist.github.com/ajpelu/194e721077ec045a2b864088908e7aca
## x is a ghlt object

table.glht <- function(x) {
    pq <- summary(x)$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    error <- attr(pq$pvalues, "error")
    pname <- switch(x$alternativ, less = paste("Pr(<", ifelse(x$df == 0, "z", "t"), ")", sep = ""), 
                    greater = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""),
                    two.sided = paste("Pr(>|",ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
    colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df == 0, "z value", "t value"), pname)
    return(mtests)
} ## end table.glht

## display N and median
## modified from : https://medium.com/@gscheithauer/how-to-add-number-of-observations-to-a-ggplot2-boxplot-b22710f7ef80
stat_box_data <- function(y, upper_limit = max(y) + 10) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y), '\n',
                    'median =', round(median(y), 1), '\n')
    )
  )
} ## end stat_box_data

## code by rnorouzian with additions by G5W (Stack Overflow)
reg <- function(x, y, col) {
    abline(lm(y ~ x), col = col)
    abline(a=0, b=1, col = "gray")
} ## reg

panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)  {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) reg(x[ok], y[ok], col.smooth)
} ## panel.lm

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
 usr <- par("usr"); on.exit(par(usr))
 par(usr = c(0, 1, 0, 1))
 r <- abs(cor(x, y))
 txt <- format(c(r, 0.123456789), digits = digits)[1]
 txt <- paste0(prefix, txt)
 text(0.5, 0.5, txt, cex = 2, font = 4)
 } ## panel.cor


## plotting multiple heatmaps
## https://support.bioconductor.org/p/87318/
grab_grob <- function(){
    grid.echo()
    grid.grab()
}

drawGridHeatmap  <- function(hm) {
    draw(hm)
    grab_grob()
}


## Create mode function
## https://www.tutorialspoint.com/r/r_mean_median_mode.htm
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
} ## getmode


## Peter Andrews arctan scale
plog <- function(x) {
    log10(1 + x^2)
}
parctan <- function(x) {
    sqrt((atan(plog(x) / plog(5))) / atan(Inf))
} ## parctan


## MAPD
mapd <- function(x, ...) {
    az <- abs(x[2:length(x)] - x[1:(length(x)-1)])
    mz <- median(az)
    sdz <- sd(az)
    cvz <- sd(az)/mean(az)
    return(c("mapd" = mz, "mapd.sd" = sdz, "mapd.cv" = cvz))
} ## mapd

## DepMap process
DepScores <- function() {
depScores <- read.csv("res/depmap/D2_combined_gene_dep_scores.csv", row.names = 1, stringsAsFactors = FALSE)
depScores <- depScores[, grep("SOFT", names(depScores))]

attr(depScores, "hgnc.symbol") <- unlist(gsub("(.*) \\(([0-9]+)\\)", "\\1", rownames(depScores)))
attr(depScores, "entrezID") <- unlist(gsub("(.*) \\(([0-9]+)\\)", "\\2", rownames(depScores)))
rownames(depScores) <- attr(depScores, "hgnc.symbol")

saveRDS(depScores, file = "res/depmap/soft_tissue_gene_dep_scores.rds")

} ## DepScores


## Process GISTIC2 output // output only grTR
## read in gistic leasions file
## provide path to gistic directory w/ all_lesions.conf_90.txt
gisticRegions <- function(gisticDir, conf = 80) {
    gr <- read.delim(file.path(gisticDir, paste0("all_lesions.conf_", conf, ".txt")), stringsAsFactors = FALSE)
    
    ## parse out genomic coordinates, and bin coordinates
    gr$chr <- gsub("chr", "", gsub("(chr[0-9XY]+):.*", "\\1", gr$Wide.Peak.Limits))
    gr$start <- gsub("(chr[0-9XY]+):([0-9]+)-([0-9]+)\\(.*", "\\2", gr$Wide.Peak.Limits)
    gr$end <- gsub("(chr[0-9XY]+):([0-9]+)-([0-9]+)\\(.*", "\\3", gr$Wide.Peak.Limits)
    gr$bin.start <- gsub("(chr[0-9XY]+):([0-9]+)-([0-9]+)\\(probes ([0-9]+):([0-9]+)\\) *", "\\4", gr$Wide.Peak.Limits)
    gr$bin.end <- gsub("(chr[0-9XY]+):([0-9]+)-([0-9]+)\\(probes ([0-9]+):([0-9]+)\\) *", "\\5", gr$Wide.Peak.Limits)
    
    ## get threshold and copy number data
    grCN <- gr[grep("Actual Copy Change Given", gr$Amplitude.Threshold),]
    grTR <- gr[grep("Actual Copy Change Given", gr$Amplitude.Threshold, invert = TRUE),]
    return(grTR)
} ## gisticRegions


## overlap gistic regions with gene / bin index
gisticGR <- function(grCN) {
    require(GenomicRanges)
    ## set up vector of amplified regions
    amp <- grCN[grep("Amp", grCN$Unique.Name), c("chr", "start", "end")]
## set up vector of deleted regions
    del <- grCN[grep("Del", grCN$Unique.Name), c("chr", "start", "end")]
    
    gr <- GRanges(
        seqnames = c(amp$chr, del$chr),
        ranges = IRanges(start = as.numeric(c(amp$start, del$start)), end = as.numeric(c(amp$end, del$end))),
        alteration.type = c(rep("Amplification", nrow(amp)), rep("Deletion", nrow(del))))
    
    geneIndex <- GRanges(seqnames = gene.index$seqnames,
                         ranges = IRanges(start = gene.index$start,

                                          end = gene.index$end),
                         strand = gene.index$strand,
                         ensembl_gene_id = gene.index$ensembl_gene_id,
                         hgnc.symbol = gene.index$hgnc.symbol,
                         gene.type = gene.index$gene.type,
                         bind.id = gene.index$bin.id)
    
    overlaps <- findOverlaps(geneIndex, gr)
    
    gistic.genes <- data.frame(gr[overlaps@to,], geneIndex[overlaps@from,])
    
    names(gistic.genes) <- c("seqnames", "chrom", "start", "end", "strand",
                             "alteration.type",
                             "chr", "gene.start", "gene.end", "width", "gene.strand",
                             "ensembl_gene_id", "hgnc.symbol", "gene.type", "bin.id")
    gistic.genes
}


## annotate with oncoKB
## give gistic genes, and oncoKB object
oncoKB_annotate <- function(gistic.genes, oncokb) {
    gistic.genes$oncokb <- "not.oncokb"
    gistic.genes$oncokb[gistic.genes$hgnc.symbol %in% oncokb$Hugo.Symbol] <- "oncokb"
    gistic.genes$oncokb[gistic.genes$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Oncogene]] <- "oncogene"
    gistic.genes$oncokb[gistic.genes$hgnc.symbol %in% oncokb$Hugo.Symbol[oncokb$Is.Tumor.Supressor.Gene]] <- "tsg"
    gistic.genes
}


## annotate essential genes w/DepMap

depmap_annotate <- function(gistic.genes, depMeans) {
## set boundaries at -0.8 for essential genes
    gistic.genes$depMeans <- depMeans[gistic.genes$hgnc.symbol]
    gistic.genes$essential <- "non-essential"
    gistic.genes$essential[gistic.genes$depMeans < -0.5] <- "essential"

    gistic.genes
    
} ## depmap annotate



#' estimate FGA
estimate_fga <- function(x, base.ploidy = 2) {
    x <- round(x)
    out <- sum(x != base.ploidy) / length(x)
    return(out)
}

