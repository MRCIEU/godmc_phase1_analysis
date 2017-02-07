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


# MR Base database
ao <- available_outcomes()
a <- subset(ao,
	priority == 1 &
	category %in% c("Disease", "Risk factor") &
	population %in% c("European", "Mixed") & 
	! id %in% c(981, 1026, 1080, 
		1089, 1090, 1031, 1085, 27, 1083, 276, 278, 1059, 1060, 991, 1074, 9, 10, 12, 980, 1013, 294, 982, 965, 966, 985, 986, 280, 286, 1025, 118, 811, 281, 283, 832, 1076, 1077, 1078, 1071, 1072, 967, 815, 32, 968, 755, 837) &
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
temp <- strsplit(mrbase$exposure, split="\\|") %>% sapply(function(x) x[1]) %>% gsub(" $", "", .)
mrbase$trait <- temp


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
g <- subset(gwas_catalog, pval < 5e-8)

g <- group_by(g, Phenotype_simple) %>%
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



# GWAS catalog
urate <- subset(gwas_catalog, Phenotype == "Urate levels (mg/dl increase)" & grepl("Kottgen", Author))
crp <- subset(gwas_catalog, grepl("C-reactive protein", Phenotype, ignore.case=TRUE) & grepl("Dehghan", Author))

gc <- format_data(rbind(urate, crp))
gc$data_source.exposure <- "GWAS catalog"

# Is there an easy way to make a sane list from GWAS catalog?

g <- subset(gwas_catalog, pval < 5e-8) %>%
	group_by(Phenotype_simple, PubmedID) %>%
	summarise(count=n()) %>%
	arrange(desc(count)) %>%
	filter(!duplicated(Phenotype_simple))

# Not really...

# eqtls

eq <- subset(gtex_eqtl, tissue == "Whole Blood" & pval < 5e-10)
eq <- format_gtex_eqtl(eq)

# proteomic qtls
pr <- format_proteomic_qtls(proteomic_qtls)

# metabolic qtls
met <- format_metab_qtls(metab_qtls)


gcc <- bind_rows(gc, eq, pr, met, mrbase)

save(gcc, file="../data/gwas_snps.rdata")


# Control SNPs - random sample of HM3 genic SNPs

hm3_genic <- scan("../data/genicsnps.txt.gz", what="character")
hm3_genic <- sample(hm3_genic, 1000, replace=FALSE)

table(hm3_genic %in% gcc$SNP)

write.table(hm3_genic, file = "../data/hm3_genic_control.txt", row=F, col=F, qu=F)



# Neanderthal SNPs
# Taken from http://science.sciencemag.org/content/351/6274/737
# Data obtained from https://phewascatalog.org/neanderthal
# Note - they quote 1495 SNPs with likely neanderthal alleles but this data source only contains 1300 (nominally significant in PHEWAS)

neanderthal <- read.csv("../data/neanderthal-phewas-catalog.csv.gz")
neanderthal <- subset(neanderthal, !duplicated(snp), select=c(chr, bp, snp, maf_eur, maf_sas, maf_afr, maf_amr, maf_eas, genes))

save(neanderthal, file="../data/neanderthal.rdata")



## Summarise

load("../data/gwas_snps.rdata")
load("../data/neanderthal.rdata")
controls <- scan("../data/hm3_genic_control.txt", what="character")


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
	return(d)
}

test <- ucsc_get_position(unique(gcc$SNP))



dat <- rbind(
	data.frame(SNP = gcc$SNP, reason = gcc$exposure, source = gcc$data_source.exposure),
	data.frame(SNP = neanderthal$snp, reason = "Neanderthal alleles", source = "Simonti et al 2016"),
	data.frame(SNP = controls, reason = "Control SNPs", source = "Randomly selected HapMap3 SNPs from genic regions")
)

group_by(dat, reason) %>%
	dplyr::summarise(n=n(), source=first(source)) %>%
	arrange(desc(source)) %>% as.data.frame()

group_by(dat, source) %>%
	dplyr::summarise(n=n())

