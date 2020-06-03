## prep fasta files for MEDICC
## imports significant cnv regions from GISTIC2.0 and calculates
## the modal copy number for that segment in each cell

## End line is two before the first number match. I.e. the line before the first appearence, and -1 for the header
grep("Actual Copy", readLines("out_WD5816_grch37.50k.k50.nobad.varbin.short.cbs_c90/all_lesions.conf_90.txt"))[1]-2

gcalls <- read.delim("out_WD5816_grch37.50k.k50.nobad.varbin.short.cbs_c90/all_lesions.conf_90.txt", as.is = TRUE, nrow = 157)

sreg <- data.frame(Wide.Peak.Limits = gsub(" +$", "", gcalls$Wide.Penak.Limits))

sreg$chr <- gsub("^chr", "", gsub("\\:.*", "", sreg$Wide.Peak.Limits))
sreg$start <- gsub("^chr.*:([0-9]+)-([0-9]+)\\(.*", "\\1", sreg$Wide.Peak.Limits)
sreg$end <- gsub("^chr.*:([0-9]+)-([0-9]+)\\(.*", "\\2", sreg$Wide.Peak.Limits)
sreg$XM1 <- as.numeric(gsub(".*\\(probes ([0-9]+):([0-9]+)\\).*", "\\1", gcalls$Peak.Limits))
sreg$X1 <- as.numeric(gsub(".*\\(probes ([0-9]+):([0-9]+)\\).*", "\\1", sreg$Wide.Peak.Limits))
sreg$X2 <- as.numeric(gsub(".*\\(probes ([0-9]+):([0-9]+)\\).*", "\\2", sreg$Wide.Peak.Limits))
sreg$nmark <- as.numeric(sreg$X2) - as.numeric(sreg$X1)

cn50k <- read.delim("vbData/WD5816_grch37.50k.k50.nobad.varbin.data.txt")
cn50k[,4:ncol(cn50k)] <- round(cn50k[,4:ncol(cn50k)])

## gCN50k <- cbind(sreg[, c("chr", "start", "end", "nmark")],  cn50k[sreg$XM1,4:ncol(cn50k)])
## gCN50k <- gCN50k[order(gCN50k$chr, gCN50k$start, gCN50k$end), ]

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

if(!dir.exists("medicc")) dir.create("medicc")

chrom <- unique(gsub(":.*", "", rownames(medicc.cnv.major)))
cells <- names(cn50k[,4:ncol(cn50k)])

sapply(chrom, function(chr) {
    ss = grep(paste0(chr, ":"), rownames(medicc.cnv.major))
    fa = as.data.frame(apply(medicc.cnv.major[ss,], 2, paste0, collapse = ""))
    rownames(fa) <- paste(">", rownames(fa))
    write.table(fa, file = file.path("medicc", paste(chr, ".cnv.major.fa", sep = "")), sep = "\n", quote = FALSE, row.names = TRUE, col.names = FALSE)
    fa1 = as.data.frame(apply(medicc.cnv.minor[ss,], 2, paste0, collapse = ""))
    rownames(fa1) <- paste(">", rownames(fa1))
    write.table(fa1, file = file.path("medicc", paste(chr, ".cnv.minor.fa", sep = "")), sep = "\n", quote = FALSE, row.names = TRUE, col.names = FALSE)
})

desc <- data.frame(chrom, major = file.path("medicc", paste(chrom, ".cnv.major.fa", sep = "")), minor = file.path("medicc", paste(chrom, ".cnv.minor.fa", sep = "")))

write.table(desc, file = file.path("medicc", "desc.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

if(!dir.exists("sifit")) dir.create("sifit")

## sifit ternary format -- separates gains from amplifications and losses, from hz deletions
im3 <- gcalls[,cells]
write.table(im, file = file.path("sifit", "gistic.calls.im3.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

## phangorn/sifit binary format plain incidence matrix: 1,0 present, absent
im <- im3
im[im != 0 ] <- 1
write.table(im, file = file.path("sifit", "gistic.calls.im.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

