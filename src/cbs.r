## Collect arguments
args <- commandArgs(trailingOnly = TRUE)
## args <-  c("--genomePath=/home/gularter/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/varbin/", "--cell.varbin=varbin50k/WD4495_35_1844.b49.dd.50k.varbin.out.txt", "--seq.stats=varbin50k/WD4495_35_1844.b49.dd.50k.varbin.stats.txt",  "--bin.size=50k",  "--ploidy.range=1.5,2.8")
 
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --genomePath=/path/to/genome/files       - character, path to varbin genome reference files
      --cell.varbin=cell-1.varbin.50k.txt      - character, name of varbin.out file with bin
      --seq.stats=cell-1.50k.varbin.stats.txt  - character, name of varbin.stats file with bin
      --bin.size=(50k|20k|5k)                  - character, bin size for the analysis
      --ploidy.range=1.5,4                       - expected ploidy range (MUST BE DIVISIBLE BY 0.05)
      --help                                   - print this text
 
      Example:
      Rscript cbs.r --genomePath=/home/gularter/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/varbin/ --cell.varbin=cell-1.varbin.50k.txt --seq.stats=cell-1.varbin.50k.stats.txt --bin.size=50k --ploidy.range=1.5,4 \n\n")
 
  q(save="no")
}

## TODO:  add min.mult and max.mult to arguments

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
argsL

## Arg1 default
if(! file.exists(argsL$genomePath)) {
    stop("genomePath ", argsL$genomePath, " not found")
    q(save = "no")
}

## Arg2 default
if(! file.exists(argsL$cell.varbin)) {
    stop("cell file ", argsL$cell.varbin, " not found")
    q(save = "no")
}

## Arg3 default
if(! argsL$bin.size %in% c("50k", "20k", "5k")) {
    stop("resources for bin.size ", argsL$bin.size, " are not available")
    q(save = "no")
}

## Arg4 default
min.ploidy <- min(strsplit(argsL$ploidy.range, split = ",")[[1]])
max.ploidy <- max(strsplit(argsL$ploidy.range, split = ",")[[1]])

## Arg5 default
if(! file.exists(argsL$seq.stats)) {
    stop("Sequencing statistics file ", argsL$seq.stats, " not found")
    q(save = "no")
}



## library
library("DNAcopy")

source("./src/cbsLib.R")

genome.path <- argsL$genomePath
cellID <- gsub("\\.[0-9]{1,2}k\\.varbin\\.out\\.txt", "", argsL$cell.varbin)


## lapply(cells, function(cellID) {
cbs.segment01(indir = ".", outdir = "./", 
              varbin.gc = file.path(genome.path, paste("varbin.gc.content.", argsL$bin.size, ".bowtie.k50.grch37.txt", sep = "")),
              varbin.data = argsL$cell.varbin, stat.data = argsL$seq.stats,
              sample.name = cellID, alt.sample.name = "",
              alpha = 0.05, nperm = 1000, undo.SD = 1.0, min.width = 5,
              mult.min = min.ploidy, mult.max = max.ploidy, bin.size = argsL$bin.size )
## })


##bad.bins = file.path(genome.path, paste("hg19.", bin.size,".k50.bad.bins.txt", sep = "")),
