library(dplyr)

load("../data/chromatin_qtls_1e-7.rdata")

table(dat$p.value < 1e-10)

filter(dat, p.value < 1e-10) %>% summary()
filter(dat, p.value < 1e-10) %>% filter(!duplicated(snpid)) %>% dim

chromatin <- subset(dat, p.value < 1e-10)
chromatin <- data.frame(id=chromatin$snpid, exposure=chromatin$peakid, pval.exposure=chromatin$p.value, stringsAsFactors=FALSE)
chromatin$data_source.exposure <- "Chromatin QTLs"

save(chromatin, file="../data/snps_chromatin.rdata")
