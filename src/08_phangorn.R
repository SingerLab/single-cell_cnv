## running phangorn
## use conda activate single-cell-cnv
## R must be the miniconda3/envs/single-cell-cnv/bin/R
R.home()

## libraries
library(ape)
library(phangorn)

## data 
file <- "phangorn/WD5816.phangorn.aa.fa"
dat = read.phyDat(file, type = "AA", format = "fasta")

## model is F81 b/c we are not using organism based aa frequency models
dm <- dist.ml(dat, model = "F81")

tree = NJ(dm)
treeUPGMA <- upgma(dm)

## testing to other models
mt <- modelTest(dat, model=c("F81", "JC69", "Blosum62"), multicore = TRUE, mc.cores = 4)

## (mt <- modelTest(dat, model="all", multicore=TRUE, mc.cores = 4))

save.image("mt.fjb.rda")

fitStart = eval(get(mt$Model[which.min(mt$BIC)], env), env)

fitNJ = pml(tree, dat, model="JTT", k=4, inv=.2)

fit = optim.pml(fitNJ, rearrangement = "stochastic",
                optInv=TRUE, optGamma=TRUE)

fit

bs = bootstrap.pml(fit, bs=1000, optNni=TRUE, multicore=TRUE, mc.cores = 4)




## init tutorial -- for DNA sequences
## fdir <- system.file("extdata/trees", package = "phangorn")
## primates <- read.phyDat(file.path(fdir, "primates.dna"), format = "interleaved")
## dm  <- dist.ml(primates)
## treeUPGMA  <- upgma(dm)
## treeNJ  <- NJ(dm)
## pdf("primate.tree.pdf")
## layout(matrix(c(1,2), 2, 1), height=c(1,2))
## par(mar = c(0,0,2,0)+ 0.1)
## plot(treeUPGMA, main="UPGMA")
## plot(treeNJ, "unrooted", main="NJ")
## dev.off()
