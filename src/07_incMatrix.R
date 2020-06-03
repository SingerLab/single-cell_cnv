## prep fasta files for MEDICC
## imports significant cnv regions from GISTIC2.0 and calculates
## the modal copy number for that segment in each cell

## End line is two before the first number match. I.e. the line before the first appearence, and -1 for the header
grep("Actual Copy", readLines("out_WD5816_grch37.50k.k50.nobad.varbin.short.cbs_c90/all_lesions.conf_90.txt"))[1]-2

gcalls <- read.delim("out_WD5816_grch37.50k.k50.nobad.varbin.short.cbs_c90/all_lesions.conf_90.txt", as.is = TRUE, nrow = 157)

sreg <- data.frame(Wide.Peak.Limits = gsub(" +$", "", gcalls$Wide.Peak.Limits))
sreg$chr <- gsub("^chr", "", gsub("\\:.*", "", sreg$Wide.Peak.Limits))
sreg$start <- gsub("^chr.*:([0-9]+)-([0-9]+)\\(.*", "\\1", sreg$Wide.Peak.Limits)
sreg$end <- gsub("^chr.*:([0-9]+)-([0-9]+)\\(.*", "\\2", sreg$Wide.Peak.Limits)
sreg$XM1 <- as.numeric(gsub(".*\\(probes ([0-9]+):([0-9]+)\\).*", "\\1", gcalls$Peak.Limits))
sreg$X1 <- as.numeric(gsub(".*\\(probes ([0-9]+):([0-9]+)\\).*", "\\1", sreg$Wide.Peak.Limits))
sreg$X2 <- as.numeric(gsub(".*\\(probes ([0-9]+):([0-9]+)\\).*", "\\2", sreg$Wide.Peak.Limits))
sreg$nmark <- as.numeric(sreg$X2) - as.numeric(sreg$X1)

cn50k <- read.delim("vbData/WD5816_grch37.50k.k50.nobad.varbin.data.txt")
cn50k[,4:ncol(cn50k)] <- round(cn50k[,4:ncol(cn50k)])

cells <- names(cn50k[,4:ncol(cn50k)])

## MEDICC input data
medicc.cnv.major <- cn50k[sreg$XM1,4:ncol(cn50k)]
rownames(medicc.cnv.major) <- gsub("\\(probes.*", "", gcalls$Peak.Limits)
any(is.na(medicc.cnv.major))

medicc.cnv.minor <- data.frame(matrix(1, nrow = nrow(medicc.cnv.major), ncol = ncol(medicc.cnv.major), dimnames = dimnames(medicc.cnv.major)))
any(is.na(medicc.cnv.minor))
sum(medicc.cnv.major <= 1)

medicc.cnv.minor[medicc.cnv.major <= 1] <- 0
sum(medicc.cnv.major <= 1)
sum(medicc.cnv.minor == 0)

medicc.cnv.major[medicc.cnv.major >= 9] <- 9
medicc.cnv.major <- medicc.cnv.major - medicc.cnv.minor

any(medicc.cnv.major > 9)

chrom <- unique(gsub(":.*", "", rownames(medicc.cnv.major)))

if(!dir.exists("medicc")) dir.create("medicc")

sapply(chrom, function(chr) {
    ss = grep(paste0(chr, ":"), rownames(medicc.cnv.major))
    fa = as.data.frame(apply(medicc.cnv.major[ss,], 2, paste0, collapse = ""))
    rownames(fa) <- paste(">", rownames(fa))
    write.table(fa, file = file.path("medicc", paste("WD5816.", chr, ".cnv.major.fa", sep = "")), sep = "\n", quote = FALSE, row.names = TRUE, col.names = FALSE)
    fa1 = as.data.frame(apply(medicc.cnv.minor[ss,], 2, paste0, collapse = ""))
    rownames(fa1) <- paste(">", rownames(fa1))
    write.table(fa1, file = file.path("medicc", paste("WD5816.", chr, ".cnv.minor.fa", sep = "")), sep = "\n", quote = FALSE, row.names = TRUE, col.names = FALSE)
})

desc <- data.frame(chrom, major = file.path("medicc", paste("WD5816.", chrom, ".cnv.major.fa", sep = "")), minor = file.path("medicc", paste(chrom, ".cnv.minor.fa", sep = "")))

write.table(desc, file = file.path("medicc", "desc.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

if(!dir.exists("sifit")) dir.create("sifit")

## sifit ternary format -- separates gains from amplifications and losses, from hz deletions
im3 <- gcalls[,cells]
rownames(im3) <- sreg$XM1
im3 <- im3[order(as.numeric(rownames(im3))), ]
write.table(im3, file = file.path("sifit", "WD5816.im3.txt"), sep = " ", col.names = FALSE, row.names = TRUE, quote = FALSE)

## phangorn/sifit binary format plain incidence matrix: 1,0 present, absent
im <- im3
im[im != 0 ] <- 1
write.table(im, file = file.path("sifit", "WD5816.im.txt"), sep = " ", col.names = FALSE, row.names = TRUE, quote = FALSE)
write(cells, file = file.path("sifit", "cellID.txt"))


## phangorn AA code change
phangorn.cnv <- cn50k[sreg$XM1,4:ncol(cn50k)]
rownames(phangorn.cnv) <- sreg$XM1
phangorn.cnv <- phangorn.cnv[order(as.numeric(rownames(phangorn.cnv))), ]

## asign a code to each value covering the complete copy number spectrum
phangorn.aa <- phangorn.cnv
phangorn.aa[phangorn.cnv == 0] <- "A"
phangorn.aa[phangorn.cnv == 1] <- "C"
phangorn.aa[phangorn.cnv == 2] <- "D"
phangorn.aa[phangorn.cnv == 3] <- "E"
phangorn.aa[phangorn.cnv == 4] <- "F"
phangorn.aa[phangorn.cnv == 5] <- "G"
phangorn.aa[phangorn.cnv == 6] <- "H"
phangorn.aa[phangorn.cnv == 7] <- "I"
phangorn.aa[phangorn.cnv == 8] <- "K"
phangorn.aa[phangorn.cnv == 9] <- "L"
phangorn.aa[phangorn.cnv >= 10 & phangorn.cnv < 20] <- "M"
phangorn.aa[phangorn.cnv >= 20 & phangorn.cnv < 30] <- "N"
phangorn.aa[phangorn.cnv >= 30 & phangorn.cnv < 40] <- "P"
phangorn.aa[phangorn.cnv >= 40 & phangorn.cnv < 50] <- "Q"
phangorn.aa[phangorn.cnv >= 50 & phangorn.cnv < 60] <- "R"
phangorn.aa[phangorn.cnv >= 60 & phangorn.cnv < 70] <- "S"
phangorn.aa[phangorn.cnv >= 70 & phangorn.cnv < 80] <- "T"
phangorn.aa[phangorn.cnv >= 80 & phangorn.cnv < 90] <- "V"
phangorn.aa[phangorn.cnv >= 90 & phangorn.cnv < 100] <- "W"
phangorn.aa[phangorn.cnv >= 100] <- "Y"

## paste as single sequence
phangorn.seq <- as.data.frame(apply(phangorn.aa, 2, paste0, collapse = ""))

## create directory
if(! dir.exists("phangorn")) dir.create("phangorn")

write.table(phangorn.seq, file = file.path("phangorn", "WD5816.phangorn.aa.txt"), quote = FALSE, row.names = TRUE, col.names = FALSE, sep = " ")
