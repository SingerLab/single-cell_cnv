## reading flow core

library(flowCore)
library(ggcyto)

fcs.dir <- dir("tmp/FCS_representative_hist")

fcs.dir <- grep("INX", fcs.dir, value = TRUE, invert = TRUE)

frames <- lapply(dir(fcs.dir, full.names=TRUE), read.FCS, transform = FALSE, alter.names = TRUE)
( fs <- as(frames, "flowSet") )

names(frames) <- sapply(frames, keyword, "SAMPLE ID")
( fs <- as(frames, "flowSet") )

sampleNames(fs)


autoplot(x, "FL1.H", "FL2.H")
autoplot(x, "FL1.H") + theme_bw()

fs <- read.flowSet(path = system.file("extdata",
                                      package = "flowCore"),
                   pattern = "\\.")

autoplot(fs, "FL1-H", "FL2-H")

gs <- GatingSet(fs)
