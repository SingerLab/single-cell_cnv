## plotting read counts

asplit <- read.table("asplit/total.read_counts.txt", header = FALSE)
bsplit <- read.table("bsplit/total.read_counts.txt", header = FALSE)
btsplit <- read.table("btsplit/total.read_counts.txt", header = FALSE)


asplit$V1 <- gsub("asplit/(.*).fq.gz", "\\1", asplit$V1)
bsplit$V1 <- gsub("bsplit/(.*).fq.gz", "\\1", bsplit$V1)
btsplit$V1 <- gsub("(.*).fq.gz", "\\1", btsplit$V1)



msplit <- merge(btsplit, bsplit, by = "V1", all = TRUE)
## msplit <- merge(msplit, asplit, by = "V1", all = TRUE)

names(msplit) <- c("btsplit", "bsplit")
