library(dplyr)


## Summarise

load("../data/snps_gwas.rdata")
load("../data/snps_metabolomic.rdata")
load("../data/snps_proteomic.rdata")
load("../data/snps_expression.rdata")
load("../data/snps_chromatin.rdata")
load("../data/snps_neanderthal.rdata")
load("../data/SOMA_pQTL.RData")
load("../data/ispc_eqtl.rdata")
w<-read.table("../data/Westra_transsnps_fdr0.05.txt",he=F)


dat <- rbind(
	data.frame(SNP = gwas$id, reason = gwas$exposure, source = gwas$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = metabolomic$id, reason = metabolomic$exposure, source = metabolomic$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = proteomic$id, reason = proteomic$exposure, source = proteomic$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = expression$id, reason = expression$exposure, source = expression$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = w$V1, reason = w$V1, source = "Westra transSNPs FDR<0.05", stringsAsFactors=FALSE),
	data.frame(SNP = chromatin$id, reason = chromatin$exposure, source = chromatin$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = p$id, reason = p$SOMAmerID, source = "proteomic_qtls", stringsAsFactors=FALSE),
	data.frame(SNP = r$SNP, reason = r$gene_id, source = "ispc_eqtl", stringsAsFactors=FALSE),
	data.frame(SNP = neanderthal$id, reason = NA, source = "Neanderthal alleles, Simonti et al 2016", stringsAsFactors=FALSE)


)

# group_by(dat, reason) %>%
# 	dplyr::summarise(n=n(), source=first(source)) %>%
# 	arrange(desc(source)) %>% as.data.frame()

sum1 <- group_by(dat, source) %>%
	filter(!duplicated(SNP)) %>%
	dplyr::summarise(nsnp=n())
sum1 <- rbind(sum1, data.frame(source="Total unique SNPs", nsnp=length(unique(dat$SNP))))

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
