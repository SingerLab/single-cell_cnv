rm(list = ls())
## run facets
R.home()

## libraries
library(facets)
source("workflows/myLib.R")

normal.sample <- "WD5816" 
tumor.sample <- "WD5816" 
inDir <- "facets/snp_pileup"


## snpMat <- file.path(inDir, paste(normal.sample, "_", tumor.sample, ".snpmat.gz", sep = "")) )
snpMat <- list.files(inDir, pattern = ".snpmat.gz", full.names = TRUE)

rcmat <- readSnpMatrix(snpMat)

cval1 <- 30
cval2 <- 400

cat("pre-cval: ", cval1 , "; proc-cval: ", cval2, "\n")
xx <- preProcSample(rcmat, het.thresh = 0.25, gbuild="hg19", cval = cval1, unmatched = TRUE)
oo <- procSample(xx, cval = cval2)
cat("proc-diplod LogR: ", oo$dipLogR, "\n")

fit <- emcncf(oo)
head(fit$cncf)
cat("tumor purity: ", fit$purity, "\n")
cat("esti'd ploidy: ", fit$ploidy, "\n")
cat("fit-diplod LogR: ", fit$dipLogR, "\n")

if(! dir.exists(paths = "facets/plots"))  dir.create(path = "facets/plots", recursive = TRUE)
if(! dir.exists(paths = "facets/cncf"))  dir.create(path = "facets/cncf", recursive = TRUE)


pdf(file.path("facets", "plots", paste(tumor.sample, "_", normal.sample, ".ploidy_", format(fit$ploidy, digits = 2), ".cnlr.pdf", sep = "")))
plotSample(x = oo, emfit = fit, sname = tumor.sample)
dev.off()

pdf(file.path("facets", "plots", paste(tumor.sample, "_", normal.sample, ".ploidy_", format(fit$ploidy, digits = 2), ".cncf.pdf", sep = "")), width = 298/25.4, height =100/25.4)
plotCNCF(x = oo, emfit = fit, main = tumor.sample)
dev.off()

write.table(fit$cncf, file = file.path("facets", "cncf", paste(tumor.sample, "_", normal.sample, ".ploidy_", format(fit$ploidy, digits = 2),".cncf.txt", sep = "")), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

fitSeg <- data.frame(ID = rep(tumor.sample, nrow(fit$cncf)))
fitSeg <- cbind(fitSeg, fit$cncf[, c("chrom", "start", "end", "num.mark", "cnlr.median")])

write.table(fit$cncf, file = file.path("facets", "cncf", paste(tumor.sample, "_", normal.sample, ".ploidy_", format(fit$ploidy, digits = 2),".cncf.seg", sep = "")), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
