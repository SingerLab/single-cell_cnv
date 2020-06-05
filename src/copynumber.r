# Collect arguments
args <- commandArgs(trailingOnly = TRUE)
 
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --nobad.wg filename     - character, name of whole genome hg19.50k.k50.nobad.varbin.data file
      --nobad.short filename  - character, name of whole genome hg19.50k.k50.nobad.varbin.short.txt
      --help                  - print this text
 
      Example:
      Rscript copynumber.r --nobad.wg=sample1.hg19.50k.k50.nobad.varbin.data.txt --nobad.short=sample1.hg19.50k.k50.nobad.varbin.short.txt \n\n")
 
  q(save="no")
}
 
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
argsL


## Arg1 default
if(! file.exists(argsL$nobad.wg)) {
    stop("`nobad.wg` file not found")
    q(save = "no")
}

## Arg2 default
if(! file.exists(argsL$nobad.short)) {
    stop("`nobad.short` file not found")
    q(save = "no")
}

source("cbsLib.R")

cellID <- gsub("(.*).hg19.*.txt", "\\1", argsL$nobad.wg)

df <- read.table(argsL$nobad.wg, header = TRUE)
dfs <- read.table(argsL$nobad.short, header = TRUE)

starts <- c()
ends <- c()
prevEnd <- 0
len <- nrow(dfs)
for (j in 1:len) {
	thisStart = prevEnd + 1
	thisEnd = thisStart + dfs$num.mark[j] - 1
	starts <- c(starts, thisStart)
	ends <- c(ends, thisEnd)
	prevEnd = thisEnd
}

amat <- matrix(data = 0, nrow = 1500000, ncol=1)
counter <- 1
for (j in 1:(len-1)) {
	for (k in (j+1):len) {
		N <- round((starts[j] - ends[j] + 1) * (starts[k] - ends[k] + 1)/1000)
		D <- abs(2^dfs$seg.mean[j] - 2^dfs$seg.mean[k])
		cat(N, "\t")
		if (N > 0) {
			amat[(counter:(counter+N-1)), 1] <- rep.int(D, N)
			counter <- counter+N
		}
	}
}
a3 <- amat[(1:counter),1]
a3.95 <- sort(a3)[round(.95*counter)]
a3d <- density(a3[which(a3 < a3.95)], n = 1000)
cn0 <- a3d$x[which(peaks(as.vector(a3d$y), span  =59))][1]
cn1 <- a3d$x[which(peaks(as.vector(a3d$y), span = 59))][2]

df$cn.ratio <- df$lowratio / cn1
df$cn.seg <- df$seg.mean.LOWESS / cn1
df$copy.number <- round(df$cn.seg)

write.table(df, sep = "\t", file = paste(cellID, ".hg19.50k.k50.varbin.data.copynumber.txt", sep = ""), quote = FALSE, row.names = FALSE)

png(paste(cellID, ".wg.cn.density.png", sep = ""), height = 148, width = 298, units = "mm", res = 350)
par(mar = c(5.1,4.1,4.1,4.1))
plot(a3d, main = paste(cellID, "seg.mean difference density"))
dev.off()

png(paste(cellID, ".wg.cn.png", sep=""), height = 148, width = 298, units = "mm", res = 350)
par(mar = c(5.1,4.1,4.1,4.1))
plot(x = df$abspos, y = df$cn.ratio, main=paste("CN profile: ", cellID, sep=""), xlab="Bin", ylab="Ratio", col="#CCCCCC")
lines(x = df$abspos, y = df$cn.ratio, col="#CCCCCC")
points(x = df$abspos, y = df$cn.seg, col="#0000DD")
lines(x = df$abspos, y = df$cn.seg, col="#0000DD")
points(x  = df$abspos, y = df$copy.number, col="#DD0000")
lines(x = df$abspos, y = df$copy.number, col="#DD0000")
dev.off()

for (a in 1:24) {
	png(paste(cellID, ".chr", a, ".cn.png", sep=""), height = 148, width = 298, units = "mm", res = 350)
	par(mar = c(5.1,4.1,4.1,4.1))
	plot(df$cn.ratio[df$chrom == a], main=paste(cellID, " chr", a, sep = ""), xlab = "Bin", ylab = "Ratio", col = "#CCCCCC")
	lines(df$cn.ratio[df$chrom == a], col = "#CCCCCC")
	points(df$cn.seg[df$chrom == a], col = "#0000DD")
	lines(df$cn.seg[df$chrom == a], col = "#0000DD")
	points(df$copy.number[df$chrom == a], col = "#DD0000")
	lines(df$copy.number[df$chrom == a], col = "#DD0000")
	dev.off()
}

