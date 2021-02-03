## Analysis with QDNAseq
## which R is being used ! must be in single-cell-cnv environment
R.home()
## Collect arguments
args <- commandArgs(trailingOnly = TRUE)
## args <- c("--sample.name=WD5816.bulk", "--input.dir=bwa_map/", "--output.dir=qdnaseq_bwa.bulk.map10", "--bin.size=50kb", "aligner=bwa")

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
      --bin.size=(50kb|100kb|500kb)                  - character, bin size in kb for the analysis
      --aligner=bowtie                               - character, bin size in kb for the analysis
      --help                                         - print this text
 
      Example:
      Rscript 05_qdnaseq.R --sample.name=cell-1 --input.dir=bowtie_map/ --output.dir=qdnaseq/ --bin.size=500k --aligner=bowtie")
 
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


## cluster for parallel comp
##options(mc.cores=8L)
##cl <- parallel::makeForkCluster(nnodes = getOption("mc.cores", 8L))
##future::plan("cluster", cluster = cl)

## Download bins or read bins if pre-downloaded
## bins <- getBinAnnotations(binSize = 50)
bins <- readRDS(paste("res/qdnaseq/qdnaseq.bins.", bin.size, ".rds", sep = ""))

## number of reads per bin -- use CACHE -- saves time later
## dev: bamList <- list.files("bowtie_map/", pattern = "dd.bam$", full.names = TRUE)[1:10]  
## dev: readCounts <- binReadCounts(bins, bamfiles = bamList , cache = TRUE)

saveRDS(readCounts, file = file.path(outDir, paste(sample.name, ".readCounts.", aligner, ".", bin.size, ".rds", sep = "")))
## readCounts <- readRDS(file.path(outDir, paste(sample.name, ".readCounts.", aligner, ".", bin.size, ".rds", sep = "")))

## flitering, normalizing and smoothing
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE, mappability = 10, chromosomes = "", )
readCountsFiltered <- estimateCorrection(readCountsFiltered)
( copyNumbers <- correctBins(readCountsFiltered) )
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

## segmentation
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt", min.width = 5)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
## yielded NA    Segmenting: WD5816_DD_1_1041.b99 (178 of 895) ..

pdf(file.path(outDir, paste(sample.name, ".copy_numbers_segmented.", aligner, ".", bin.size, ".pdf", sep = ""), width = 210/25.4, height = 148/25.4)
plot(copyNumbersSegmented, pointcol = "gray", ylim = c(-5, 5))
dev.off()

copyNumbersCalled <- callBins(copyNumbersSegmented, method = "CGHcall")
saveRDS(copyNumbersCalled, file = file.path(outDir, paste(sample.name, ".copy_number_calls.", aligner, ".", bin.size, ".rds", sep = "")))

pdf(file.path(outDir, paste(sample.name, ".copy_numbers_called.", aligner, ".", bin.size, ".pdf", sep = "")), width = 210/25.4, height = 148/25.4)
plot(copyNumbersCalled, pointcol = "gray", ylim = c(-5, 5))
dev.off()

## ( cgh <- makeCgh(copyNumbersCalled) )

save.image(file.path(outDir, paste(sample.name, aligner, bin.size, "qdnaseq.out.rda", sep = ".")))

#exportBins(copyNumbersCalled, format = "vcf")
exportBins(copyNumbersCalled, file = file.path(outDir, paste(sample.name, ".single-cell.", aligner, ".", bin.size, "_calls.txt", sep = "")), type = "calls", logTransform = FALSE)  ## OK

exportBins(copyNumbersCalled, file = file.path(outDir, paste(sample.name, ".single-cell.", aligner, ".", bin.size, "_integer_copynumber.txt", sep = "")), type = "copynumber", logTransform = FALSE)  ## seems ok

exportBins(copyNumbersCalled, file = file.path(outDir, paste(sample.name, ".single-cell.", aligner, ".", bin.size, "_log2_copynumber.txt", sep = "")), type = "copynumber", logTransform = TRUE)


tryCatch(exportBins(copyNumbersCalled, format = "seg", includeZero = TRUE))
## source("src/exportSEG.w0.R")
## tryCatch(exportSEGw0(copyNumbersCalled, includeZero = TRUE))

system(paste("cat *seg | awk '/SAMPLE_NAME/ NR > 1 {next} {print $0}' | sed -e 's/.sorted//' >", file.path(outDir, paste(sample.name, ".single-cell.", aligner, ".", bin.size,  ".igv.seg", sep = ""))))
list.files(path = outDir, pattern = "seg$")
system(paste("rm ", sample.name, "*.seg", sep = ""))


## QDNAseq:::betterCall()
## pdf(file.path(outDir, paste(sample.name, ".read_counts.", aligner, ".", bin.size, ".pdf", sep = "")), width = 210/25.4, height = 148/25.4)
## nplot(readCounts, logTransform = FALSE, ylim = c(-50, 200))
## highlightFilters(readCounts, logTransform = FALSE, residual = TRUE, blacklist = TRUE)
## dev.off()
## pdf(file.path(outDir, paste(sample.name, ".copy_numbers_smooth.", aligner, ".", bin.size, ".pdf", sep = "")), width = 210/25.4, height = 148/25.4)
## plot(copyNumbersSmooth)
## dev.off()
## exportBins(copyNumbersSmooth, file = file.path(outDir, paste(sample.name, ".", aligner, ".", bin.size, "_smothCN_single-cell.txt", sep = "")))
## exportBins(copyNumbersSmooth, file = file.path(outDir, paste(sample.name, ".", aligner, ".", bin.size, "_smothCN_single-cell.igv", sep = "")), format = "igv")
## exportBins(copyNumbersSmooth, file = file.path(outDir, paste(sample.name, ".", aligner, ".", bin.size, "_smothCN_single-cell.bed", sep = "")), format = "bed")

