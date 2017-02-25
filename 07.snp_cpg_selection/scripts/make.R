library(stringr)
library(dplyr)
library(tidyr)
library(meffil)


# CpGs

load("../data/cpgs_variable.rdata")

feat <- meffil.get.features()

cpglist <- unique(cpgdat$cpg)
a <- feat[feat$name %in% cpglist, ]
a$chromosome <- gsub("X", "23", a$chromosome)
a$chromosome <- gsub("Y", "24", a$chromosome)
a$chr <- as.numeric(gsub("chr", "", a$chromosome))
b <- arrange(a, chr, position)

cpglist <- b$name

chunk <- function(x, n) split(x, sort(rank(x) %% n))
c <- chunk(cpglist, 100)

for(i in 1:100)
{
	write.table(c[[i]], file=paste0("../lists/cpglist_", i, ".txt"), row=F, col=F, qu=F)
}


# SNPs

load("../data/snplist.rdata")

# exclude controls
dat <- subset(dat, ! source %in% c("Randomly selected HapMap3 SNPs from genic regions", "SNPs with evidence for selection, high global Fst and MAF in Europeans"))
dat <- subset(dat, !duplicated(SNP))

a <- str_split_fixed(dat$SNP, ":", 3) %>% as_data_frame()
a$V1 <- gsub("X", "23", a$V1)
a$V1 <- gsub("Y", "24", a$V1)
a$V1 <- gsub("chr17_ctg5_hap1", "chr17", a$V1)
a$V1 <- gsub("chr6_cox_hap2", "chr6", a$V1)
a$chr <- as.numeric(gsub("chr", "", a$V1))
a$pos <- as.numeric(a$V2)
a$id <- dat$SNP

a <- a[order(a$chr, a$pos), ]
a$id <- gsub("X", "23", a$id)

write.table(a$id, file="../lists/snplist.txt", row=F, col=F, qu=F)
