## Jude Kendall functions for single-cell copy number analysis
require("DNAcopy")

lowess.gc <- function(jtkx, jtky) {
    jtklow <- lowess(jtkx, log(jtky), f = 0.05)
    jtkz <- approx(jtklow$x, jtklow$y, jtkx)
    return(exp(log(jtky) - jtkz$y))
} ## end lowess.gc

remove.segment <- function(rsShort, rsSegnum, ratioData, sd.undo) {

    appendLeft <- TRUE
    checkSdundo <- FALSE

    if (rsSegnum == 1) {
        appendLeft <- FALSE
    } else {
	if (rsSegnum == nrow(rsShort)) {
            appendLeft <- TRUE
	} else {
            rightIndex <- rsSegnum + 1
            leftIndex <- rsSegnum - 1
            
            if (rsShort[rightIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
                appendLeft <- TRUE
            } else {
		if (rsShort[leftIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
                    appendLeft <- FALSE
		} else {
                    if (abs(rsShort[leftIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"]) < abs(rsShort[rightIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"])) {
			appendLeft <- TRUE
			checkSdundo <- TRUE
                    } else {
			appendLeft <- FALSE
			checkSdundo <- TRUE
                    }}}
	}}
    
    appendIndex <- 99999999
    if (appendLeft) {
        appendIndex <- rsSegnum - 1
    } else {
        appendIndex <- rsSegnum + 1
    }
    
    tempShort <- rsShort
    newLocStart <- -1
    newLocEnd <- -1
    if (appendLeft) {
        tempShort[appendIndex, "loc.end"] <- tempShort[rsSegnum, "loc.end"]
        tempShort[appendIndex, "seg.end"] <- tempShort[rsSegnum, "seg.end"]
    } else {
        tempShort[appendIndex, "loc.start"] <- tempShort[rsSegnum, "loc.start"]
        tempShort[appendIndex, "seg.start"] <- tempShort[rsSegnum, "seg.start"]
    }
    
    tempShort[appendIndex, "num.mark"] <- tempShort[appendIndex, "num.mark"] + tempShort[rsSegnum, "num.mark"]
    tempShort[appendIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[appendIndex, "seg.start"]:tempShort[appendIndex, "seg.end"]], base = 2))
    
    cat("append", tempShort[appendIndex, "chrom"], tempShort[appendIndex, "loc.start"], tempShort[appendIndex, "loc.end"], tempShort[appendIndex, "num.mark"], tempShort[appendIndex, "seg.mean"], tempShort[appendIndex, "seg.start"], tempShort[appendIndex, "seg.end"], "\n")
    
    tempShort <- tempShort[-rsSegnum, ]
    tempShort$segnum <- seq(1:nrow(tempShort))
    
    if (checkSdundo) {
        thisSd <- -1
        if (appendLeft) {
            leftIndex <- appendIndex
            rightIndex <- appendIndex + 1
        } else {
            leftIndex <- appendIndex - 2
            rightIndex <- appendIndex - 1
        }
                                        #thisSd <- sd(ratioData[tempShort$seg.start[leftIndex]:tempShort$seg.start[rightIndex], "lowratio"])
        thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)
        
        if (abs(tempShort$seg.mean[leftIndex] - tempShort$seg.mean[rightIndex]) < (sd.undo * thisSd) ) {

            cat("left", tempShort[leftIndex, "chrom"], tempShort[leftIndex, "loc.start"], tempShort[leftIndex, "loc.end"], tempShort[leftIndex, "num.mark"], tempShort[leftIndex, "seg.mean"], tempShort[leftIndex, "seg.start"], tempShort[leftIndex, "seg.end"], "\n")
            cat("right", tempShort[rightIndex, "chrom"], tempShort[rightIndex, "loc.start"], tempShort[rightIndex, "loc.end"], tempShort[rightIndex, "num.mark"], tempShort[rightIndex, "seg.mean"], tempShort[rightIndex, "seg.start"], tempShort[rightIndex, "seg.end"], "\n")

            ## remove breakpoint
            tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
            tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
            tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
            tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base = 2))
            tempShort <- tempShort[-rightIndex, ]
            tempShort$segnum <- seq(1:nrow(tempShort))
        }
    }
    
    return(tempShort)
} ## end remove.segment


sdundo.all <- function (sdShort, ratioData, sd.undo) {

    tempShort <- sdShort
    thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)
    
    while ( TRUE ) {
        
        chrom <- tempShort$chrom
        chrom.shift <- c(tempShort$chrom[-1], tempShort$chrom[1])
	
        breakpoints <- which(chrom == chrom.shift)
        cat("sdundo.all intrachrom breakpoints", length(breakpoints), "\n")
        
        if (length(breakpoints) < 1) {
            break
        }
        
        breakpoints.shift <- breakpoints + 1
        
        undo.breakpoints <- breakpoints[which(abs(tempShort$seg.mean[breakpoints] - tempShort$seg.mean[breakpoints.shift]) < thisSd * sd.undo)]

        cat("sdundo.all undo breakpoints", length(undo.breakpoints), "\n")

        if (length(undo.breakpoints) < 1) {
            break
        }
        
        undo.breakpoints.shift <- undo.breakpoints + 1
        
        undo.df <- tempShort[undo.breakpoints, ]
        undo.df$seg.mean.diff <- abs(tempShort$seg.mean[undo.breakpoints] - tempShort$seg.mean[undo.breakpoints.shift])

        min.index <- which.min(undo.df$seg.mean.diff)
        
        leftIndex <- undo.df$segnum[min.index]
        rightIndex <- leftIndex + 1

        cat("sdundo.all left", tempShort[leftIndex, "chrom"], tempShort[leftIndex, "loc.start"], tempShort[leftIndex, "loc.end"], tempShort[leftIndex, "num.mark"], tempShort[leftIndex, "seg.mean"], tempShort[leftIndex, "seg.start"], tempShort[leftIndex, "seg.end"], "\n")
        cat("sdundo.all right", tempShort[rightIndex, "chrom"], tempShort[rightIndex, "loc.start"], tempShort[rightIndex, "loc.end"], tempShort[rightIndex, "num.mark"], tempShort[rightIndex, "seg.mean"], tempShort[rightIndex, "seg.start"], tempShort[rightIndex, "seg.end"], "\n")

        tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
        tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
        tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
        tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base = 2))
        tempShort <- tempShort[-rightIndex, ]
        tempShort$segnum <- seq(1:nrow(tempShort))

    }
    
    return(tempShort)
    
} ## end sdundo.all

numeric.chrom <- function(chrom) {
    nch <- gsub("chr", "", chrom)
    nch[as.character(nch) == "X"] <- 23
    nch[as.character(nch) == "Y"] <- 24
    nch[as.character(nch) %in% c("M", "MT")] <- 25
    as.numeric(nch)
} ## end numeric.chr

mapd <- function(x, ...) {
    az <- abs(x[2:length(x)] - x[1:(length(x)-1)])
    mz <- median(az)
    sdz <- sd(az)
    cvz <- sd(az)/mean(az)
    return(c("mapd" = mz, "mapd.sd" = sdz, "mapd.cv" = cvz))
} ## end mapd

cbs.segment01 <- function(indir, outdir, varbin.gc, varbin.data, stat.data, sample.name, alt.sample.name, alpha, nperm, undo.SD, min.width, mult.min, mult.max, bin.size) {

    gc <- read.table(varbin.gc, header = TRUE)

    ## get vector of chromosomes
    chrom.numeric <- numeric.chrom(gc$bin.chrom) ## see function above
    
    thisRatio <- read.table(paste(indir, varbin.data, sep = "/"), header = FALSE)
    names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
    thisRatio$chrom <- numeric.chrom(thisRatio$chrom)
    thisRatio <- thisRatio[order(thisRatio$chrom, thisRatio$chrompos),]
    thisRatio$abspos <- cumsum(as.numeric(thisRatio$chrompos))
    a <- thisRatio$bincount + 1
    thisRatio$ratio <- a / mean(a)
    thisRatio$gc.content <- gc$gc.content
    thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)
    ##thisRatio$log2lowratio <- log(thisRatio$lowratio, base = 2) 

    mapd.qc <- mapd(thisRatio$lowratio)
    
    set.seed(81)
    CNA.object <- CNA(log(thisRatio$lowratio, base = 2), thisRatio$chrom, thisRatio$chrompos, data.type = "logratio", sampleid = sample.name) 
    smoothed.CNA.object <- smooth.CNA(CNA.object) 
    segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha = alpha, nperm = nperm, undo.splits = "sdundo", undo.SD = undo.SD, min.width = 2)
    thisShort <- segment.smoothed.CNA.object[[2]]

    m <- matrix(data = 0, nrow = nrow(thisRatio), ncol = 1)	
    prevEnd <- 0
    for (i in 1:nrow(thisShort)) {
        thisStart <- prevEnd + 1
        thisEnd <- prevEnd + thisShort$num.mark[i]
        m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
        prevEnd = thisEnd
    }
    cbs.long <- m[, 1]

##### NEW STUFF also check min.width = 2 above
    
    write.table(thisShort, sep = "\t", file = file.path(outdir, paste(sample.name, ".grch37.", bin.size, ".k50.varbin.short.cbs.seg", sep = "")), quote = FALSE, row.names = FALSE) 

    workShort <- thisShort
    workShort$segnum <- 0
    workShort$seg.start <- 0
    workShort$seg.end <- 0
    prevEnd <- 0
    for (i in 1:nrow(thisShort)) {
        thisStart <- prevEnd + 1
        thisEnd <- prevEnd + thisShort$num.mark[i]
        workShort$seg.start[i] <- thisStart
        workShort$seg.end[i] <- thisEnd
        workShort$segnum[i] <- i
        prevEnd = thisEnd
    }

    discardSegments <- TRUE
    while (discardSegments) {
        orderShort <- workShort[order(workShort$num.mark, abs(workShort$seg.mean)), ]
        if (orderShort[1, "num.mark"] < min.width) {
            workShort <- remove.segment(workShort, orderShort[1, "segnum"], thisRatio, undo.SD)
        } else {
            discardSegments <- FALSE
        }
    }

    workShort <- sdundo.all(workShort, thisRatio, undo.SD)
    thisShort <- workShort 
    
##### END NEW STUFF
    
    m <- matrix(data = 0, nrow = nrow(thisRatio), ncol = 1)	
    prevEnd <- 0
    for (i in 1:nrow(thisShort)) {
        thisStart <- prevEnd + 1
        thisEnd <- prevEnd + thisShort$num.mark[i]
        m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
        prevEnd = thisEnd
    }
    thisRatio$seg.mean.LOWESS <- m[, 1]

    thisGrid <- seq(mult.min, mult.max, by = 0.05)
    thisOuter <- thisRatio$seg.mean.LOWESS %o% thisGrid
    thisOuterRound <- round(thisOuter)
    thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
    thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
    thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
    thisError <- min(thisOuterColsums)
    thisShredded <- length(which(thisRatio$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

    thisRatio$ratio.quantal <- thisRatio$lowratio * thisMultiplier
    thisRatio$seg.quantal <- thisRatio$seg.mean.LOWESS * thisMultiplier
    
    thisRatio$cbs.seg <- cbs.long
    thisRatio$cbs.seg.quantal <- cbs.long * thisMultiplier

    thisQuantalStats <- data.frame("ploidy" = thisMultiplier, "error" = thisError, "shredded" = thisShredded)
    
    chr <- thisRatio$chrom
    chr.shift <- c(chr[-1], chr[length(chr)])
    vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
    hlines <- c(0.5, 1.0, 1.5, 2.0)
    chr.text <- c(1:22, "X", "Y", "MT")
    vlines.shift <- c(vlines[-1], 4*10^9)
    chr.at <- vlines + (vlines.shift - vlines) / 2
    x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
    x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

    ## add stats -- ploidy is now added to the plot
    ## add stats -- total number of reads / duplicates
    bam.stats <- read.delim(file.path(indir, stat.data), as.is = TRUE)
    bam.stats$GenomeCoverage <- ((bam.stats$ReadsKept * 52) / 3098825702)
    write.table(bam.stats, file = file.path(indir, stat.data), sep = "\t", quote = FALSE, row.names = FALSE)
    
    tryCatch({
        ## pdf(file.path(outdir, paste(sample.name, ".", bin.size, ".wg.pdf", sep = "")), height = 149/25.4, width = 298/25.4)
        png(file.path(outdir, paste(sample.name, ".", bin.size, ".wg.png", sep = "")), height = 149, width = 298, units = "mm", res = 350)
        plot(x = thisRatio$abspos, y = thisRatio$lowratio,
             log = "y", main = paste(sample.name, alt.sample.name),
             xaxt = "n", xlab = "Genome Position Gb", ylab = "Ratio", col = "#CCCCCC", cex = 0.4)
        lines(x = thisRatio$abspos, y = thisRatio$lowratio, col = "#CCCCCC", lwd = 0.8)
        points(x = thisRatio$abspos, y = thisRatio$seg.mean.LOWESS, col = "#D40000", cex = 0.8)
        lines(x = thisRatio$abspos, y = thisRatio$seg.mean.LOWESS, col="#D40000")
        abline(h = hlines, lty = 3)
        abline(v = vlines, lty = 3)
        axis(1, at = chr.at, labels = chr.text, cex = 0.8, tick = FALSE)
        mtext(text = paste("Ploidy =", thisMultiplier, "; ",
                           "Total Reads =", bam.stats$TotalReads, "; ",
                           "Reads Kept =", bam.stats$ReadsKept, "; ",
                           "Median Bin Count =", bam.stats$MedianBinCount, "; ",
                           "MAPD =", round(mapd.qc["mapd"],2), "+/-", round(mapd.qc["mapd.sd"], 2), "; ",
                           "Avg. Genome Coverage =", format(bam.stats$GenomeCoverage, digits = 2)), adj = 0)
        dev.off()
    })
    hlines <- c(1, 2, 3, 4, 5, 6)

    tryCatch({
        ## pdf(file.path(outdir, paste(sample.name, ".", bin.size, ".wg.quantal.pdf", sep = "")), height = 149/25.4, width = 298/25.4)
        png(file.path(outdir, paste(sample.name, ".", bin.size, ".wg.quantal.png", sep = "")), height = 149, width = 298, units = "mm", res = 350)
        plot(x = thisRatio$abspos, y = thisRatio$ratio.quantal,
             log = "y", main = paste(sample.name, alt.sample.name), xaxt = "n",
             xlab = "Genome Position (Chr)", ylab = "Ratio Quantal", col = "#CCCCCC", cex = 0.4)
        lines(x = thisRatio$abspos, y = thisRatio$ratio.quantal, col = "#CCCCCC", lwd = 0.8)
        points(x = thisRatio$abspos, y = thisRatio$seg.quantal, col = "#D40000", cex = 0.8)
        lines(x = thisRatio$abspos, y = thisRatio$seg.quantal, col = "#D40000")
        abline(h = hlines, lty = 3)
        abline(v = vlines, lty = 3)
        axis(1, at = chr.at, labels = chr.text, cex = 0.8, tick = FALSE)
        mtext(text = paste("Ploidy =", thisMultiplier, "; ",
                           "Total Reads =", bam.stats$TotalReads, "; ",
                           "Reads Kept =", bam.stats$ReadsKept, "; ",
                           "Median Bin Count =", bam.stats$MedianBinCount, "; ",
                           "MAPD =", round(mapd.qc[["mapd"]], 2), "+/-", round(mapd.qc[["mapd.sd"]], 2), "; ",
                           "Avg. Genome Coverage =", format(bam.stats$GenomeCoverage, digits = 2)), adj = 0)
        dev.off()
    })
    
    write.table(thisQuantalStats, sep = "\t", file = file.path(outdir, paste(sample.name, ".grch37.", bin.size, ".k50.varbin.quantal.stats.txt", sep = "")), quote = FALSE, row.names = FALSE) 
    write.table(thisRatio, sep = "\t", file = file.path(outdir, paste(sample.name, ".grch37.", bin.size, ".k50.varbin.data.txt", sep="")), quote = FALSE, row.names = FALSE) 
    write.table(thisShort, sep = "\t", file = file.path(outdir, paste(sample.name, ".grch37.", bin.size, ".k50.varbin.short.txt", sep="")), quote = FALSE, row.names = FALSE) 
    
    ## re-analysis removing bad bins
    bad <- read.table("/home/gularter/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/varbin/grch37.50k.k50.bad.bins.txt",
                      header = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
    
    thisRatioNobig <- thisRatio[-bad[, 1], ]
    
    set.seed(81)
    CNA.object <- CNA(log(thisRatioNobig$lowratio, base = 2), thisRatioNobig$chrom, thisRatioNobig$chrompos, data.type = "logratio", sampleid = sample.name) 
    smoothed.CNA.object <- smooth.CNA(CNA.object) 
    segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha = alpha, nperm = nperm, undo.splits = "sdundo", undo.SD = undo.SD, min.width = 2) 
    thisShort <- segment.smoothed.CNA.object[[2]]
    
    m <- matrix(data = 0, nrow = nrow(thisRatioNobig), ncol = 1)	
    prevEnd <- 0
    for (i in 1:nrow(thisShort)) {
        thisStart <- prevEnd + 1
        thisEnd <- prevEnd + thisShort$num.mark[i]
        m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
        prevEnd = thisEnd
    }
    cbs.long.nobad <- m[, 1]
    
##### NEW STUFF also check min.width=2 above

##%     write.table(thisShort, sep = "\t", file = file.path(outdir, paste(sample.name, ".grch37.", bin.size, ".k50.nobad.varbin.short.cbs.seg", sep = "")), quote = FALSE, row.names = FALSE) 
    
    workShort <- thisShort
    workShort$segnum <- 0
    workShort$seg.start <- 0
    workShort$seg.end <- 0
    prevEnd <- 0
    for (i in 1:nrow(thisShort)) {
        thisStart <- prevEnd + 1
        thisEnd <- prevEnd + thisShort$num.mark[i]
        workShort$seg.start[i] <- thisStart
        workShort$seg.end[i] <- thisEnd
        workShort$segnum[i] <- i
        prevEnd = thisEnd
    }
    
    discardSegments <- TRUE
    while (discardSegments) {
        orderShort <- workShort[order(workShort$num.mark, abs(workShort$seg.mean)), ]
        if (orderShort[1, "num.mark"] < min.width) {
            workShort <- remove.segment(workShort, orderShort[1, "segnum"], thisRatioNobig, undo.SD)
        } else {
            discardSegments <- FALSE
        }
    }
    
    workShort <- sdundo.all(workShort, thisRatioNobig, undo.SD)
    thisShort <- workShort
    
##### END NEW STUFF
    
    
    m <- matrix(data = 0, nrow = nrow(thisRatioNobig), ncol = 1)	
    prevEnd <- 0
    for (i in 1:nrow(thisShort)) {
        thisStart <- prevEnd + 1
        thisEnd <- prevEnd + thisShort$num.mark[i]
        m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
        prevEnd = thisEnd
    }
    
    thisRatioNobig$seg.mean.LOWESS <- m[, 1]
    
    thisGrid <- seq(mult.min, mult.max, by = 0.05)
    thisOuter <- thisRatioNobig$seg.mean.LOWESS %o% thisGrid
    thisOuterRound <- round(thisOuter)
    thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
    thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
    thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
    thisError <- min(thisOuterColsums)
    thisShredded <- length(which(thisRatioNobig$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

    thisRatioNobig$ratio.quantal <- thisRatioNobig$lowratio * thisMultiplier
    thisRatioNobig$seg.quantal <- thisRatioNobig$seg.mean.LOWESS * thisMultiplier

    thisRatioNobig$cbs.seg <- cbs.long.nobad
    thisRatioNobig$cbs.seg.quantal <- cbs.long.nobad * thisMultiplier
    
    thisQuantalStatsNobad <- data.frame("ploidy" = thisMultiplier, "error" = thisError, "shredded" = thisShredded)
    
    chr <- thisRatioNobig$chrom
    chr.shift <- c(chr[-1], chr[length(chr)])
    ## Use the vlines abspos positions from above. Because these start at
    ## the third bin on the acrocentric chromosomes the vlines end up to
    ## the right of the centromere rather than the left which is wrong.
    ## vlines <- c(1, thisRatioNobig$abspos[which(chr != chr.shift) + 1], thisRatioNobig$abspos[nrow(thisRatioNobig)])
    hlines <- c(0.5, 1.0, 1.5, 2.0)
    chr.text <- c(1:22, "X", "Y", "MT")
    vlines.shift <- c(vlines[-1], 4*10^9)
    chr.at <- vlines + (vlines.shift - vlines) / 2
    x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
    x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

##%     tryCatch({
##%         ## pdf(file.path(outdir, paste(sample.name, ".", bin.size, ".wg.nobad.pdf", sep = "")), height = 149/25.4, width = 298/25.4)
##%         png(file.path(outdir, paste(sample.name, ".", bin.size, ".wg.nobad.png", sep = "")), height = 149, width = 298, units = "mm", res = 350)
##%         plot(x = thisRatioNobig$abspos, y = thisRatioNobig$lowratio,
##%           log = "y", main = paste(sample.name, alt.sample.name),
##%              xaxt = "n", xlab = "Genome Position Gb", ylab = "Ratio", col = "#CCCCCC", cex = 0.4)
##%         lines(x = thisRatioNobig$abspos, y = thisRatioNobig$lowratio, col = "#CCCCCC", lwd = 0.8)
##%         points(x = thisRatioNobig$abspos, y = thisRatioNobig$seg.mean.LOWESS, col = "#D40000", cex = 0.8)
##%         lines(x = thisRatioNobig$abspos, y = thisRatioNobig$seg.mean.LOWESS, col = "#D40000")
##%         axis(1, at = chr.at, labels = chr.text, cex = 0.8, tick = FALSE)
##%         abline(h = hlines, lty = 3)
##%         abline(v = vlines, lty = 3)
##%         mtext(text = paste("Ploidy =", thisMultiplier, "; ",
##%                            "Total Reads =", bam.stats$TotalReads, "; ",
##%                            "Reads Kept =", bam.stats$ReadsKept, "; ",
##%                            "Median Bin Count =", bam.stats$MedianBinCount, "; ",
##%                            "MAPD =", round(mapd.qc["mapd"],2), "+/-", round(mapd.qc["mapd.sd"], 2), "; ",
##%                            "Avg. Genome Coverage =", format(bam.stats$GenomeCoverage, digits = 2)), adj = 0)
##%         dev.off()
##%     })
    
    hlines <- c(1, 2, 3, 4, 5, 6)

##%     tryCatch({
##%         pdf(file.path(outdir, paste(sample.name, ".", bin.size, ".wg.nobad.quantal.pdf", sep = "")), height = 149/25.4, width = 298/25.4)
##%          png(file.path(outdir, paste(sample.name, ".", bin.size, ".wg.nobad.quantal.png", sep = "")), height = 149, width = 298, units = "mm", res = 350) 
##%         plot(x = thisRatioNobig$abspos, y = thisRatioNobig$ratio.quantal,
##%              log = "y", main = paste(sample.name, alt.sample.name),
##%              xaxt = "n", xlab = "Genome Position Gb", ylab = "Ratio Quantal", col = "#CCCCCC", cex = 0.4)
##%         axis(1, at = chr.at, labels = chr.text, cex = 0.8, tick = FALSE)
##%         lines(x = thisRatioNobig$abspos, y = thisRatioNobig$ratio.quantal, col = "#CCCCCC", lwd = 0.8)
##%         points(x = thisRatioNobig$abspos, y = thisRatioNobig$seg.quantal, col = "#D40000", cex = 0.8)
##%         lines(x = thisRatioNobig$abspos, y = thisRatioNobig$seg.quantal, col = "#D40000")
##%         abline(h = hlines, lty = 3)
##%         abline(v = vlines, lty = 3)
##%         mtext(text = paste("Ploidy =", thisMultiplier, "; ",
##%                            "Total Reads =", bam.stats$TotalReads, "; ",
##%                            "Reads Kept =", bam.stats$ReadsKept, "; ",
##%                            "Median Bin Count =", bam.stats$MedianBinCount, "; ",
##%                            "MAPD =", round(mapd.qc["mapd"],2), "+/-", round(mapd.qc["mapd.sd"], 2), "; ",
##%                            "Avg. Genome Coverage =", format(bam.stats$GenomeCoverage, digits = 2)), adj = 0)
##%         dev.off()
##%     })
    
##%     write.table(thisQuantalStatsNobad, sep = "\t", file = file.path(outdir, paste(sample.name, ".grch37.", bin.size, ".k50.nobad.varbin.quantal.stats.txt", sep = "")), quote = FALSE, row.names = FALSE) 
##%     write.table(thisRatioNobig, sep = "\t", file = file.path(outdir, paste(sample.name, ".grch37.", bin.size, ".k50.nobad.varbin.data.txt", sep = "")), quote = FALSE, row.names = FALSE) 
##%     write.table(thisShort, sep = "\t", file = file.path(outdir, paste(sample.name, ".grch37.", bin.size, ".k50.nobad.varbin.short.txt", sep = "")), quote = FALSE, row.names = FALSE) 
    
    
} ## end cbs.segment01



peaks <- function(series, span = 3, ties.method = "first") {
### This peaks function from Petr Pikal https://stat.ethz.ch/pipermail/r-help/2007-February/125241.html
    if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
    z <- embed(series, span)
    s <- span%/%2
    v <- max.col(z, ties.method = ties.method) == 1 + s
    pad <- rep(FALSE, s)
    result <- c(pad, v, pad)
    result
} ## end peaks

