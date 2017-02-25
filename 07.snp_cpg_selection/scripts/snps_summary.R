library(dplyr)




## Summarise

load("../data/snps_gwas.rdata")
load("../data/snps_metabolomic.rdata")
load("../data/snps_proteomic.rdata")
load("../data/snps_expression.rdata")
load("../data/snps_neanderthal.rdata")
load("../data/snps_selection.rdata")
controls <- read.table("../data/snps_control.txt", header=FALSE)

dat <- rbind(
	data.frame(SNP = gwas$id, reason = gwas$exposure, source = gwas$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = metabolomic$id, reason = metabolomic$exposure, source = metabolomic$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = proteomic$id, reason = proteomic$exposure, source = proteomic$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = expression$id, reason = expression$exposure, source = expression$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = neanderthal$id, reason = NA, source = "Neanderthal alleles, Simonti et al 2016", stringsAsFactors=FALSE),
	data.frame(SNP = as.character(controls$V2), reason = NA, source = "Randomly selected HapMap3 SNPs from genic regions", stringsAsFactors=FALSE),
	data.frame(SNP = as.character(selection$SNP), reason = NA, source = "SNPs with evidence for selection, high global Fst and MAF in Europeans", stringsAsFactors=FALSE)
)

# group_by(dat, reason) %>%
# 	dplyr::summarise(n=n(), source=first(source)) %>%
# 	arrange(desc(source)) %>% as.data.frame()

sum1 <- group_by(dat, source) %>%
	dplyr::summarise(nsnp=n())

sum2 <- filter(dat, source == "mrbase") %>%
	group_by(reason) %>%
	dplyr::summarise(nsnp=n())

sum3 <- filter(dat, source == "GWAS catalog") %>%
	group_by(reason) %>%
	dplyr::summarise(nsnp=n())

save(dat, file="../data/snplist.rdata")

write.csv(sum1, "../data/sum1.csv")
write.csv(sum2, "../data/sum2.csv")
write.csv(sum3, "../data/sum3.csv")

write.table(unique(dat$SNP), file="../data/snplist.txt", row=F, col=F, qu=F)
