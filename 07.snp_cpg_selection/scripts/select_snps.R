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
	! id %in% c(981, 1026, 1080)
)
a$chunk <- rep(1:ceiling(nrow(a)/10), each=10)[1:nrow(a)]

l <- list()
for(i in 21:23)
{
	l[[i]] <- extract_instruments(a$id[a$chunk == i])
}

mrbase <- bind_rows(l)


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