exportSEGw0 <- function(obj, fnames=NULL , includeZero = FALSE) {

    require(Biobase) || stop("Biobase not installed")
    
    calls <- assayDataElement(obj, "calls")
    segments <- QDNAseq:::log2adhoc(assayDataElement(obj, "segmented"))

    fd <- fData(obj)
    pd <- pData(obj)
  
    if (is.null(fnames)) 
	fnames <- pd$name

    if (length(fnames) != length(pd$name)) {
        print("Length of names is too short")
    }
 
    for (i in 1:ncol(calls)) {	
	d <- cbind(fd[,1:3],calls[,i], segments[,i])

        if(includeZero) {
            sel <- !is.na(d[,4])
        } else {
            sel <- d[,4] != 0 & !is.na(d[,4])
        }
        
            dsel <- d[sel,]
        
	rle(paste(d[sel,1], d[sel,4], sep=":")) -> rleD

	endI <- cumsum(rleD$lengths)
	posI <- c(1, endI[-length(endI)] + 1)

	chr <- dsel[posI,1]
	pos <- dsel[posI,2]
	end <- dsel[endI,3]
	score <- dsel[posI,4]
	segVal <- round(dsel[posI,5],2)
	bins <- rleD$lengths

	options(scipen=100)

	out <- cbind(fnames[i], chr, pos, end, bins, segVal)
	colnames(out) <- c("SAMPLE_NAME", "CHROMOSOME", "START", "STOP", "DATAPOINTS", "LOG2_RATIO_MEAN")

	fname <- paste(fnames[i], ".seg", sep="")

	write.table(out, fname, quote=F, sep="\t", append=FALSE, col.names=TRUE, row.names=FALSE)
    }
    
}


    
