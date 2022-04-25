## Analysis with QDNAseq
## which R is being used ! must be in single-cell-cnv environment
R.home()
## Collect arguments
args <- commandArgs(trailingOnly = TRUE)
## args <- c("--sample.name=LS8817", "--input.dir=bowtie_out/", "--output.dir=~/.Rcache", "--bin.size=50kb", "--aligner=bowtie")

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R/QDNASeq Script
 
      Arguments:
      --sample.name=cell-1                           - character, path to varbin genome reference files
      --input.dir=bowtie_map/                        - character, name of input directory with all bam files
      --output.dir=qdnaseq_bowtieXY.500kb_map10/     - character, name of outputdirectory for files to be deposited
      --bin.size=(50kb|500kb)                        - character, bin size in kb for the analysis
      --aligner=bowtie                               - character, bin size in kb for the analysis
      --help                                         - print this text
 
      Example:
      Rscript 05_qdnaseq.R --sample.name=cell-1 --input.dir=bowtie_map/ --output.dir=qdnaseq/ --bin.size=50k --aligner=bowtie")
 
  q(save="no")
}
 
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
argsL

## Arg1 default
if(! dir.exists(argsL$input.dir)) {
    stop("input directory ", argsL$input.dir, " not found")
    q(save = "no")
}

## Arg2 default -- create output directory if it doesn't exist
if(! dir.exists(argsL$output.dir)) dir.create(argsL$output.dir)


## Arg3 default
if(! argsL$bin.size %in% c("500kb", "50kb")) {
    stop("resources for bin.size", argsL$bin.size, "are not available")
    q(save = "no")
}

## libraries
library(QDNAseq)

## Run parameters
sample.name <- argsL$sample.name
inDir  <- file.path(argsL$input.dir)
outDir <- file.path(argsL$output.dir)
bin.size <- argsL$bin.size
aligner <- argsL$aligner

bins <- readRDS(paste("res/qdnaseq/qdnaseq.bins.", bin.size, ".rds", sep = ""))
readCounts <- binReadCounts(bins, path = inDir, ext = "md.bam", cache = TRUE)

## bamList <- list.files(inDir, pattern = "md.bam$")
## readCounts <- binReadCounts(bins, bamfiles = bamList , cache = TRUE)
