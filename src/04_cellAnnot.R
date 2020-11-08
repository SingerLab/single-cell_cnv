## build a simple cell annotation file
R.home() == "/work/singer/opt/miniconda3/envs/single-cell-cnv/lib/R" || stop("Wrong environment, run `conda activate single-cell-cnv`")

library(RColorBrewer)
library(wesanderson)
library(viridis)

## Run parameters
sample.name <- "WD5816"
inDir  <- file.path("vbData/")
outDir <- file.path("vbData/")
bin.size <- "5k"
aligner <- "bowtie"

hh <- as.character(read.delim(file.path(inDir, paste(sample.name, "_grch37.", bin.size, ".k50.nobad.varbin.data.txt", sep = "")), header = FALSE, as.is = TRUE)[1,])

## Cell Annotation -- best to read as ploidy estimates have ben also incorporated
cellAnnot <- data.frame(cellID = hh[4:length(hh)])
cellAnnot$subsample <- gsub("\\.b[0-9]+.*", "", unlist(sapply(strsplit(as.character(cellAnnot$cellID), split = "_"), function(i) paste(i[2], i[3], sep = "_"))))
cellAnnot$subtype <- gsub("([A-Z]{1,3})[0-9]+.*", "\\1", cellAnnot$cellID)
cellAnnot$plate <- gsub("([A-Z]{1,3}[0-9]+.*)\\.b[0-9]+.*", "\\1", cellAnnot$cellID)
cellAnnot$barcode <- gsub(".*\\.(b[0-9]+).*", "\\1", cellAnnot$cellID)

cellAnnot$plate[cellAnnot$plate == "WD5816_DD_2_1048.08" & cellAnnot$barcode %in% paste0("b", 1:96)] <- "WD5816_DD_2_1048"
cellAnnot$plate[cellAnnot$plate == "WD5816_DD_2_1048.08" & cellAnnot$barcode %in% paste0("b", 97:192)] <- "WD5816_DD_2_1108"
cellAnnot$plate[cellAnnot$plate == "WD5816_DD_3_1121.27" & cellAnnot$barcode %in% paste0("b", 1:96)] <- "WD5816_DD_3_1121"
cellAnnot$plate[cellAnnot$plate == "WD5816_DD_3_1121.27" & cellAnnot$barcode %in% paste0("b", 97:192)] <- "WD5816_DD_3_1127"


## ploidy
ploidy.5k <- read.delim(file.path(inDir, paste0(sample.name, "_grch37.", bin.size, ".k50.nobad.varbin.quantal.stats.txt")), as.is = TRUE)
ploidy.5k$cellID <- gsub("-", ".", ploidy.5k$cellID)
( names(ploidy.5k)[2:4] <- paste(names(ploidy.5k)[2:4], bin.size, sep = ".") )

cellAnnot <- merge(cellAnnot, ploidy.5k, by = "cellID")
rownames(cellAnnot) <- cellAnnot$cellID

cellAnnot$subtype <- "DDLS w/RMS"
cellAnnot$subtype[cellAnnot$subsample == "DD_2"] <- "DDLS"
## cellAnnot$subtype[cellAnnot$subsample %in% c("DD_4", "LL_2", "LL_3")] <- "UNK"

## color choices
subtypeColor <- brewer.pal(8, "Dark2")[1:length(unique(cellAnnot$subtype))]
names(subtypeColor) <- unique(cellAnnot$subtype)
cellAnnot$subtypeColor <- subtypeColor[cellAnnot$subtype]
saveRDS(subtypeColor, file = file.path(outDir, "subtypeColor.rds"))

set.seed(81)
( subsampleColor <- sample(colors(distinct = TRUE), length(unique(cellAnnot$subsample))) )
## rainbow(length(unique(cellAnnot$subsample)))
names(subsampleColor) <- unique(cellAnnot$subsample)
cellAnnot$subsampleColor <- subsampleColor[cellAnnot$subsample]
saveRDS(subsampleColor, file = file.path(outDir, "subsampleColor.rds"))

## plateColor <- brewer.pal(10, "Paired")[1:length(unique(cellAnnot$plate))]
plateColor <- rainbow(length(unique(cellAnnot$plate)))
names(plateColor) <- unique(cellAnnot$plate)
cellAnnot$plateColor <- plateColor[cellAnnot$plate]
saveRDS(plateColor, file = file.path(outDir, "plateColor.rds"))

write.table(cellAnnot, file = file.path(outDir, paste(sample.name, "_cellAnnotations.txt", sep = "")), sep = "\t", quote = TRUE, row.names = TRUE)
