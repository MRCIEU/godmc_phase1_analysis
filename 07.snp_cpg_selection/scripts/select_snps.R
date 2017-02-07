library(TwoSampleMR)
library(MRInstruments)
library(dplyr)

# longevity

# anthropometric
# - height
# - bmi
# - lung function
# - sleep

# cancer

# autoimmune / inflammatory
# - crp
# - il6
# - eczema
# - psoriasis
# - crohns
# - ulcerative colitis
# - 

# psychiatric
# - schizophrenia
# - bipolar
# - 

# personality
# - education attainment
# - neuroticism


# cardiometabolic
# - blood pressure
# - metabolic
# - chd
# - ldl
# - hdl
# - trig
# - t1d
# - t2d


ucsc_get_position <- function(snp)
{
	snp <- paste(snp, collapse="', '")
	require(RMySQL)
	message("Connecting to UCSC MySQL database")
	mydb <- dbConnect(MySQL(), user="genome", dbname="hg19", host="genome-mysql.cse.ucsc.edu")

	query <- paste0(
		"SELECT * from snp144 where name in ('", snp, "');"
	)
	# message(query)
	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)
	# dbClearResult(dbListResults(mydb)[[1]])
	dbDisconnect(mydb)

	d <- subset(d, !duplicated(name))
	d$id <- paste0(d$chrom, ":", d$chromEnd, ":SNP")
	d <- data.frame(SNP=d$name, id=d$id, stringsAsFactors=FALSE)
	return(d)
}


# MR Base database
ao <- available_outcomes()
a <- subset(ao,
	priority == 1 &
	category %in% c("Disease", "Risk factor") &
	population %in% c("European", "Mixed") & 
	! id %in% c(981, 1026, 1080
		# , 1089, 1090, 1031, 1085, 27, 1083, 276, 278, 1059, 1060, 991, 1074, 9, 10, 12, 980, 1013, 294, 982, 965, 966, 985, 986, 280, 286, 1025, 118, 811, 281, 283, 832, 1076, 1077, 1078, 1071, 1072, 967, 815, 32, 968, 755, 837
	) &
	! author %in% c("Cousminer", "Huffman JE")
)
a$chunk <- rep(1:ceiling(nrow(a)/10), each=10)[1:nrow(a)]

subset(a, select=c(trait, nsnp, sample_size, id)) %>% arrange(trait)


l <- list()
for(i in unique(a$chunk))
{
	l[[i]] <- extract_instruments(a$id[a$chunk == i])
}

mrbase <- bind_rows(l)

mrbase <- subset(gcc, data_source.exposure=="mrbase")
mrbase$trait <- strsplit(mrbase$exposure, split="\\|") %>% sapply(function(x) x[1]) %>% gsub(" $", "", .)


out <- group_by(mrbase, trait) %>%
	do({
		x <- .
		x <- arrange(x, pval.exposure)
		x <- subset(x, !duplicated(SNP))
		if(length(unique(x$exposure)) > 1)
		{
			message("clumping ", x$trait[1])
			x <- clump_data(x)
		} else {
			message("not clumping ", x$trait[1])
		}
		return(x)
	})



data(gwas_catalog)
g <- filter(gwas_catalog, pval < 5e-8) %>%
	group_by(Phenotype_simple) %>%
	do({
		x <- .
		x$pval.exposure <- x$pval
		x <- arrange(x, pval.exposure)
		x <- subset(x, !duplicated(SNP))
		x$id.exposure <- round(runif(1) * 1000000)
		if(length(unique(x$exposure)) > 1)
		{
			message("clumping ", x$Phenotype_simple[1])
			x <- clump_data(x)
		} else {
			message("not clumping ", x$Phenotype_simple[1])
		}
		return(x)
	})

gwascat <- format_data(g)
gwascat$data_source.exposure <- "GWAS catalog"


# GWAS catalog
# urate <- subset(gwas_catalog, Phenotype == "Urate levels (mg/dl increase)" & grepl("Kottgen", Author))
# crp <- subset(gwas_catalog, grepl("C-reactive protein", Phenotype, ignore.case=TRUE) & grepl("Dehghan", Author))

# gc <- format_data(rbind(urate, crp))
# gc$data_source.exposure <- "GWAS catalog"

# # Is there an easy way to make a sane list from GWAS catalog?

# g <- subset(gwas_catalog, pval < 5e-8) %>%
# 	group_by(Phenotype_simple, PubmedID) %>%
# 	summarise(count=n()) %>%
# 	arrange(desc(count)) %>%
# 	filter(!duplicated(Phenotype_simple))

# Not really...

# eqtls

eq <- subset(gtex_eqtl, tissue == "Whole Blood" & pval < 5e-10)
eq <- format_gtex_eqtl(eq)

# proteomic qtls
pr <- format_proteomic_qtls(proteomic_qtls)

# metabolic qtls
met <- format_metab_qtls(metab_qtls)


gcc <- bind_rows(gwascat, eq, pr, met, out)
snpid <- ucsc_get_position(unique(gcc$SNP))

gcc <- merge(gcc, snpid)

save(gcc, file="../data/gwas_snps.rdata")


# Control SNPs - random sample of HM3 genic SNPs

hm3_genic <- scan("../data/genicsnps.txt.gz", what="character")
hm3_genic <- sample(hm3_genic, 2000, replace=FALSE)
table(hm3_genic %in% gcc$SNP)
hm3_genic <- ucsc_get_position(hm3_genic)
hm3_genic <- hm3_genic[! hm3_genic$SNP %in% gcc$SNP, ][1:1000, ]

write.table(hm3_genic, file = "../data/hm3_genic_control.txt", row=F, col=F, qu=F)



# Neanderthal SNPs
# Taken from http://science.sciencemag.org/content/351/6274/737
# Data obtained from https://phewascatalog.org/neanderthal
# Note - they quote 1495 SNPs with likely neanderthal alleles but this data source only contains 1300 (nominally significant in PHEWAS)

neanderthal <- read.csv("../data/neanderthal-phewas-catalog.csv.gz")
neanderthal <- subset(neanderthal, !duplicated(snp), select=c(chr, bp, snp, maf_eur, maf_sas, maf_afr, maf_amr, maf_eas, genes))
nsn <- ucsc_get_position(neanderthal$snp)
neanderthal <- merge(nsn, neanderthal, by.y="snp", by.x="SNP")
save(neanderthal, file="../data/neanderthal.rdata")



## Selection SNPs

load("~/data/1000g_reference/1000GP_Phase3/vcf/pop_freq/all_fst.rdata")
selsnps <- filter(dat, EUR > 0.1, EUR < 0.9, fst > 0.5)
save(selsnps, file="../data/selection.rdata")




## Summarise

load("../data/gwas_snps.rdata")
load("../data/neanderthal.rdata")
controls <- read.table("../data/hm3_genic_control.txt", header=FALSE)

dat <- rbind(
	data.frame(SNP = gcc$id, reason = gcc$exposure, source = gcc$data_source.exposure, stringsAsFactors=FALSE),
	data.frame(SNP = neanderthal$id, reason = NA, source = "Neanderthal alleles, Simonti et al 2016", stringsAsFactors=FALSE),
	data.frame(SNP = as.character(controls$V2), reason = NA, source = "Randomly selected HapMap3 SNPs from genic regions", stringsAsFactors=FALSE),
	data.frame(SNP = as.character(selsnps$SNP), reason = NA, source = "SNPs with evidence for selection, high global Fst and MAF in Europeans", stringsAsFactors=FALSE)
)

# group_by(dat, reason) %>%
# 	dplyr::summarise(n=n(), source=first(source)) %>%
# 	arrange(desc(source)) %>% as.data.frame()

group_by(dat, source) %>%
	dplyr::summarise(n=n())

save(dat, file="../data/snplist.rdata")



### 

## Selection SNPs

test <- ucsc_get_position(unique(out$SNP))
test <- subset(test, !duplicated(name))
snp <- paste0(test$chrom, ":", test$chromEnd, ":SNP")
temp <- subset(dat, SNP %in% snp)
t.test(temp$fst[1:1000], dat$fst[sample(1:nrow(dat), 5000, replace=FALSE)])


data(aries_mqtl)
mqtl <- sample(unique(aries_mqtl$SNP), 5000, replace=FALSE)
mqtl2 <- ucsc_get_position(mqtl)
mqtl <- paste0(mqtl2$chrom, ":", mqtl2$chromEnd, ":SNP")
temp <- subset(dat, SNP %in% mqtl)
t.test(temp$fst[1:1000], dat$fst)

controls2 <- ucsc_get_position(controls)
controls <- paste0(controls2$chrom, ":", controls2$chromEnd, ":SNP")
temp <- subset(dat, SNP %in% controls)
t.test(temp$fst, dat$fst)


# It looks like mqtls, GWAS SNPs and hm3 genic control SNPs all have higher fst than average fst across the genome. 

# Just choose top fst SNPs that have maf > 0.1 in europeans

hist(selsnps$EUR)



