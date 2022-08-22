#!/opt/common/CentOS_7/R/R-4.0.0/lib64/Rscript
## matrix data and seg files from varbin 
## which R is being used ! must be in single-cell-cnv environment
## R.home() == "/opt/common/CentOS_7/R/R-4.0.0/lib64/R" || stop("Wrong environment, run `module load R`")
source("src/myLib.R")

## Collect arguments
args <- commandArgs(trailingOnly = TRUE)
## args <- c("--sample.name=WD1544", "--input.dir=varbin5k/", "--output.dir=vbData/", "--bin.size=5k", "--aligner=bowtie")

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
cell.files <- list.files(path = inDir,
                         pattern = paste(".grch37.", bin.size,
                                         ".k50.varbin.data.txt", sep = ""))

xx <- do.call(rbind, strsplit(cell.files, split = "\\."))
## reconstruct cell id
cell.list <- paste(xx[,1], xx[,2], sep = ".")
                  
## get coordinate data
chh = read.table(file.path(inDir, cell.files[1]),
                 header = TRUE)[,1:3]

cat("reading copy number varbin data data...\n")
## get seg.quantal data for each cell
vbd <- lapply(cell.files, function(cell) {
    vbseg = read.table(file.path(inDir, cell),
                       header = TRUE)
    vbseg
})
names(vbd) <- cell.list

rm(xx)

#' vbd is a list object containing the *varbin.data.txt
#' names on the list are cellID's
#' also appends 
select_vbd_column <- function(vbd, column, coord) {
    cc <- do.call(cbind, lapply(vbd, function(i) i[, column]))
    
    if(!is.null(coord)) {
        out <- data.frame(coord, cc)
    } else {
        out <- cc
    }
    return(out)
}

cat("selecting and writing quantals")
## cbs.seg.quantal is the same as seg.quantal
## chose cbs.seg.quantal colum as it's more descritive of the processing
csq.df <- select_vbd_column(vbd = vbd, column = "cbs.seg.quantal",
                                        coord = chh)
write.table(csq.df,
            file = gzfile(file.path(outDir, paste(sample.name, "_grch37.", bin.size,
                                           ".k50.varbin.data.txt.gz", sep = ""))),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
rm(csq.df)

## seg.mean.LOWESS is the same as cbs.seg;
## using seg.mean.LOWESS as it is more descriptive of column content
## get lowess ratio (lowratio) data (segmented)
cbd.df <- select_vbd_column(vbd, column = "seg.mean.LOWESS",
                            coord = chh)
write.table(cbd.df,
            file = gzfile(file.path(outDir,
                             paste(sample.name, "_grch37.", bin.size,
                                   ".k50.varbin.lowess.data.txt.gz", sep = ""))),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
rm(cbd.df)

## BIN COUNT DATA: (raw counts + gc-corrected, and unsegmented!)
## cat("reading raw bin counts...\n")
rbc.df <- select_vbd_column(vbd = vbd, column = "bincount",
                            coord = chh)
write.table(rbc.df,
            file = gzfile(file.path(outDir,
                                    paste(sample.name, "_grch37.", bin.size,
                                          ".k50.bin_counts.txt.gz", sep = ""))),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
rm(rbc.df)

## MAPD
cat("estimating MAPD...\n")
## get lowess ration (lowratio) data (not segmented)
lbd <- select_vbd_column(vbd = vbd, column = "lowratio",
                         coord = NULL)
## calculate mapd
mapd.qc <- t(apply(lbd, 2, mapd))
mapd.qc <- data.frame(cellID = rownames(mapd.qc), mapd.qc)
## write out
write.table(mapd.qc,
            file = file.path(outDir,
                             paste(sample.name, "_grch37.", bin.size,
                                   ".k50.varbin.mapd.qc.txt", sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("reading ploidy and stats data...\n")
## write out cell annotations file
q.stats.files <- list.files(inDir, pattern = "k50.varbin.quantal.stats.txt$",
                            full.names = TRUE)

xx <- do.call(rbind, strsplit(gsub(paste0("varbin", bin.size, "/"), "",
                                   q.stats.files), "\\."))
cell.list.qs <- paste(xx[,1], xx[,2], sep = ".")

quantal.stats <- do.call(rbind, lapply(q.stats.files, read.delim, nrow = 1))
quantal.stats <- data.frame(cellID = cell.list.qs, quantal.stats)

rownames(quantal.stats) <- quantal.stats$cellID

write.table(quantal.stats,
            file = file.path(outDir,
                             paste(sample.name, "_grch37.", bin.size,
                                   ".k50.varbin.quantal.ploidy.txt", sep = "")),
            sep = "\t", quote = FALSE , row.names = FALSE)
rm(quantal.stats, xx, cell.list.qs)

## get varbin stats
cat("reading varbin stats data...\n")
## write out cell annotations file
vb.stats.files <- list.files(inDir, pattern = ".varbin.stats.txt",
                             full.names = TRUE)

xx <- do.call(rbind, strsplit(gsub(paste0("varbin", bin.size, "/"), "",
                                   vb.stats.files), "\\."))
cell.list.vb <- paste(xx[,1], xx[,2], sep = ".")

varbin.stats <- do.call(rbind, lapply(vb.stats.files, read.delim, nrow = 1))
varbin.stats <- data.frame(cellID = cell.list.vb, varbin.stats)

rownames(varbin.stats) <- varbin.stats$cellID
write.table(varbin.stats, file = file.path(outDir, paste(sample.name, "_grch37.varbin.stats.txt", sep = "")),
            sep = "\t", quote = FALSE , row.names = FALSE)

## READ WRITE SEG FILES
cat("reading segment data...\n")
## write out .seg file
seg.files <- list.files(inDir, pattern = "varbin.short.cbs.seg$",
                        full.names = TRUE)

igv.seg <- do.call(rbind, lapply(seg.files, read.delim))
## names(igv.seg) <- c("ID", "chrom", "start", "end", "nmark", "ratio")

## strsplit on "\\." and paste 1 and 2 to re-build cell names
## remove varbin{5,20,50}k folder name, and use '.' regexp to remove any
## following character (either / or .)
xx <- do.call(rbind, strsplit(
                         gsub(paste0("varbin", bin.size, "."), "",
                              igv.seg$ID), "\\."))

igv.seg$ID <- paste(xx[,1], xx[,2], sep = ".")

write.table(igv.seg,
            file = file.path(outDir,
                             paste(sample.name, "_grch37.", bin.size,
                                   ".k50.varbin.short.cbs.seg", sep = "")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# extract gc.content from the first file
#% gc.content <- select_vbd_column(vbd = vbd, column = "gc.content")
## scale to geometric mean
#% TARGET_SCALE <- 2^mean(log2(colSums(rbc.df[,-c(1:3)]]+1)))
#% ## run lowess/loess normalization
#% cat("gc corrected bin counts...\n")
#% gc.norm.bin.count <- apply(rbc.df[, -c(1:3)], 2, function(rbc) {
#%     ## estimate bin weights based on gc (4 sig. digits)
#%     wts <- tapply(rbc, round(gc.content, digits = 4), sum)
#%     wts <- wts[as.character(round(gc.content, 4))]
#%     ## log read count
#%     ## TARGET_SCALE <- 1e6
#%     scl <- TARGET_SCALE / sum(rbc +1)
#%     logtn <- log2((rbc + 1) * scl)
#%     ## weighed loess regression
#%     gcb <- loess(logtn ~ gc.content, weights = wts, span = 0.3,
#%                  control = loess.control(surface = "direct"), 
#%                  degree = 2)
#%     ## scaling read bin counts
#%     res <- 2^-gcb$fitted * rbc
#%     res[is.na(res)] <- 0
#%     ## return res
#%     return(res)
#% })
#% 
#% write.table(gc.norm.bin.count,
#%             file = file.path(outDir, paste(sample.name, "_grch37.", bin.size,
#%                                            ".k50.gc_corrected.bin_counts.txt",
#%                                            sep = "")),
#%             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#% 
