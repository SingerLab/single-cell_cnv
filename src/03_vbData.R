#!/opt/common/CentOS_7/R/R-4.0.0/lib64/Rscript
## matrix data and seg files from varbin 
## which R is being used ! must be in single-cell-cnv environment
## R.home() == "/opt/common/CentOS_7/R/R-4.0.0/lib64/R" || stop("Wrong environment, run `module load R`")

source("src/myLib.R")

## Collect arguments
args <- commandArgs(trailingOnly = TRUE)
## args <- c("--sample.name=WD9048R_MR_T", "--input.dir=varbin20k/", "--output.dir=vbData/", "--bin.size=20k", "--aligner=bowtie")

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R/Varbin Tables Script
 
      Arguments:
      --sample.name=bioID                    - character, path to varbin genome reference files
      --input.dir=varbin5k                   - character, path to varbin genome reference files
      --output.dir=vbData                    - character, path to varbin genome reference files
      --bin.size=(5k,20,50k)                 - character, total bins to compile 5k, 20k, 50k
      --aligner=bowtie                       - character, bin size in kb for the analysis
      --help                                 - print this text
 
      Example:
      Rscript src/06_vbData.R --sample.name=WD5816 --input.dir=varbin50k/ --output.dir=vbData/ --bin.size=50k --aligner=bowtie
      Rscript src/06_vbData.R --sample.name=WD5816 --input.dir=varbin20k/ --output.dir=vbData/ --bin.size=20k --aligner=bowtie
      Rscript src/06_vbData.R --sample.name=WD5816 --input.dir=varbin5k/ --output.dir=vbData/ --bin.size=5k --aligner=bowtie
")
 
  q(save="no")
}


## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## sandardizing to no dashes and 'k' = thousand bins; not kilobases (kb)
argsL$input.dir <- gsub("\\/+", "", argsL$input.dir)
argsL$output.dir <- gsub("\\/+", "", argsL$output.dir)
argsL$bin.size <- gsub("kb", "k", argsL$bin.size)

argsL

## Arg1 default
if(! dir.exists(argsL$input.dir)) {
    stop("input directory ", argsL$input.dir, " not found")
    q(save = "no")
}

## Arg2 default -- create output directory if it doesn't exist
if(! dir.exists(argsL$output.dir)) {dir.create(argsL$output.dir)}


## Arg3 default
if(! argsL$bin.size %in% c("50k", "20k", "5k")) {
    stop("resources for bin.size", argsL$bin.size, "are not available")
    q(save = "no")
}


## retrieving data into matrices
## Run parameters
sample.name <- argsL$sample.name
inDir  <- file.path(argsL$input.dir)
outDir <- file.path(argsL$output.dir)
bin.size <- argsL$bin.size
aligner <- argsL$aligner

## load copy number data
cell.list <- gsub(paste0(".*(",sample.name,".*)\\.dd.*"), "\\1",
                  list.files(path = inDir,
                             pattern = paste(".dd.grch37.", bin.size, ".k50.varbin.data.txt", sep = "")))

## get coordinate data
chh = read.table(file.path(inDir, paste(cell.list[1], ".dd.grch37.", bin.size,
                                        ".k50.varbin.data.txt", sep = "")),
                 header = TRUE)[,1:3]

cat("reading copy number quantal data...\n")
## get seg.quantal data for each cell
vbd <- sapply(cell.list, function(cell) {
    vbseg = read.table(file.path(inDir, paste(cell, ".dd.grch37.", bin.size,
                                              ".k50.varbin.data.txt", sep = "")),
                       header = TRUE)[,"seg.quantal"]
    vbseg
})
## colnames(vbd) <- gsub("\\.seg\\.quantal", "", colnames(vbd))

cn <- data.frame(chh, vbd)

write.table(cn, file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                               ".k50.varbin.data.txt", sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## get cbs.seg data
cbd <- sapply(cell.list, function(cell) {
    cbseg = read.table(file.path(inDir, paste(cell, ".dd.grch37.", bin.size,
                                              ".k50.varbin.data.txt", sep = "")),
                       header = TRUE)[, "cbs.seg"]
    cbseg
})
## colnames(cbd) <- gsub("\\.cbs\\.seg", "", colnames(cbd))

cn.seg <- data.frame(chh, cbd)

write.table(cn.seg, file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                                   ".k50.varbin.cbs_seg.data.txt", sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## get lowess ratio (lowratio) data (segmented)
sbd <- sapply(cell.list, function(cell) {
    cbseg = read.table(file.path(inDir, paste(cell, ".dd.grch37.", bin.size,
                                              ".k50.varbin.data.txt", sep = "")),
                       header = TRUE)[, "seg.mean.LOWESS"]
    cbseg
})
## colnames(cbd) <- gsub("\\.cbs\\.seg", "", colnames(cbd))

seg.mean.low <- data.frame(chh, sbd)

write.table(seg.mean.low, file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                                   ".k50.varbin.lowratio.data.txt", sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


## bincount data (raw counts and unsegmented !)
cat("reading raw bin counts...\n")
bcnt <- sapply(cell.list, function(cell) {
    useg = read.table(file.path(inDir, paste(cell, ".dd.grch37.", bin.size,
                                              ".k50.varbin.data.txt", sep = "")),
                       header = TRUE)[, "bincount"]
    useg
})
## colnames(cbd) <- gsub("\\.cbs\\.seg", "", colnames(cbd))
raw.bin.count <- data.frame(chh, bcnt)

write.table(raw.bin.count,
            file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                           ".k50.bin_counts.txt", sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# extract gc.content from the first file
gc.content <- read.table(file.path(inDir, paste(cell.list[1], ".dd.grch37.",
                                                bin.size,
                                                ".k50.varbin.data.txt", sep = "")),
                         header = TRUE)[, "gc.content"]

## scale to geometric mean
TARGET_SCALE <- 2^mean(log2(colSums(bcnt+1)))

## run lowess/loess normalization
cat("gc corrected bin counts...\n")
gc.norm.bin.count <- apply(bcnt, 2, function(rbc) {
    ## estimate bin weights based on gc (4 sig. digits)
    wts <- tapply(rbc, round(gc.content, digits = 4), sum)
    wts <- wts[as.character(round(gc.content, 4))]
    ## log read count
    ## TARGET_SCALE <- 1e6
    scl <- TARGET_SCALE / sum(rbc +1)
    logtn <- log2((rbc + 1) * scl)
    ## weighed loess regression
    gcb <- loess(logtn ~ gc.content, weights = wts, span = 0.3,
                 control = loess.control(surface = "direct"), 
                 degree = 2)
    ## scaling read bin counts
    res <- 2^-gcb$fitted * rbc
    ## return res
    return(res)
})

write.table(gc.norm.bin.count,
            file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                           ".k50.gc_corrected.bin_counts.txt",
                                           sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## MAPD
cat("estimating MAPD...\n")
## get lowess ration (lowratio) data (not segmented)
lbd <- sapply(cell.list, function(cell) {
    cbseg = read.table(file.path(inDir, paste(cell, ".dd.grch37.", bin.size,
                                              ".k50.varbin.data.txt", sep = "")),
                       header = TRUE)[, "lowratio"]
    cbseg
})

mapd.qc <- t(apply(lbd, 2, mapd))
mapd.qc <- data.frame(cellID = rownames(mapd.qc), mapd.qc)

write.table(mapd.qc, file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                                    ".k50.varbin.mapd.qc.txt", sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


## READ WRITE SEG FILES
cat("reading segment data...\n")
## write out .seg file
seg.files <- list.files(inDir, pattern = "nobad.varbin.short.cbs.seg", full.names = TRUE)
igv.seg <- do.call(rbind, lapply(seg.files, read.delim))
## names(igv.seg) <- c("ID", "chrom", "start", "end", "nmark", "ratio")

igv.seg$ID <- gsub(paste0(".*(",sample.name,".*)\\.dd.*"), "\\1", igv.seg$ID)

write.table(igv.seg, file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                                    ".k50.varbin.short.cbs.seg", sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


cat("reading ploidy and stats data...\n")
## write out cell annotations file
q.stats.files <- list.files(inDir, pattern = "k50.varbin.quantal.stats.txt", full.names = TRUE)
quantal.stats <- do.call(rbind, lapply(q.stats.files, read.delim))
quantal.stats <- data.frame(cellID = gsub(paste0(".*(",sample.name,".*)\\.dd.*"), "\\1", q.stats.files), quantal.stats)
rownames(quantal.stats) <- quantal.stats$cellID

write.table(quantal.stats,
            file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                           ".k50.varbin.quantal.ploidy.txt", sep = "")),
            sep = "\t", quote = FALSE , row.names = FALSE)


## get varbin stats
cat("reading varbin stats data...\n")
## write out cell annotations file
vb.stats.files <- list.files(inDir, pattern = ".varbin.stats.txt", full.names = TRUE)
varbin.stats <- do.call(rbind, lapply(vb.stats.files, read.delim))
varbin.stats <- data.frame(cellID = gsub(paste0(".*(",sample.name,".*)\\.dd.*"), "\\1", vb.stats.files), varbin.stats)
rownames(varbin.stats) <- varbin.stats$cellID
write.table(varbin.stats, file = file.path(outDir, paste(sample.name, "_grch37.varbin.stats.txt", sep = "")),
            sep = "\t", quote = FALSE , row.names = FALSE)
