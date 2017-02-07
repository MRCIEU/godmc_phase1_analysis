library(dplyr)
library(tidyr)
library(meffil)

cmd <- "ls /panfs/panasas01/shared-godmc/sftp/GoDMC/*/*/results/01/methylation_summary.RData > filelist.txt"
system(cmd)

filelist <- scan("filelist.txt", what="character")

info <- gsub("/panfs/panasas01/shared-godmc/sftp/GoDMC/", "", filelist) %>% gsub("/results/01/methylation_summary.RData", "", .) %>%
	strsplit(., split="/") %>%
	do.call(rbind, args=.) %>%
	as.data.frame() %>%
	transmute(analyst=V1, cohort=V2, filename=filelist)


l <- list()
for(i in 1:length(filelist))
{
	message(i)
	load(filelist[i])
	meth_summary$cohort <- info$cohort[i]
	l[[i]] <- meth_summary
}

dat <- bind_rows(l)

means <- spread(subset(dat, select=c(cpg, cohort, mean)), key=cohort, value=mean)
sds <- spread(subset(dat, select=c(cpg, cohort, sd)), key=cohort, value=sd)

cormeans <- cor(means[,-c(1)], use="pair")
corsds <- cor(sds[,-c(1)], use="pair")


top_sds <- group_by(dat, cohort) %>%
	do({
		x <- .
		thresh <- quantile(x$sd, 0.7)
		return(subset(x, sd >= thresh, select=c(cpg, sd)))
	})

tab_cpg <- table(top_sds$cpg)

cpglist <- names(tab_cpg)[tab_cpg >= 23]
length(cpglist)

## 

cpgdat <- data.frame(cpg = cpglist, source = "20k most variable sites amongst all 23 GoDMC phase 1 cohorts", stringsAsFactors=FALSE)

##


feat <- meffil.get.features()

# Randomly sampled from each chromosome

feats <- subset(feat, name %in% cpglist)

a <- table(feat$chromosome) / nrow(feat)
b <- table(feats$chromosome) / nrow(feats)

plot(as.numeric(a)[1:22], as.numeric(b), xlab='Proportion of CpGs on chromosome', ylab="Proportion of selected CpGs on chromosome")


# Selected has more from open sea and shores, less from islands

# Most tissue / cancer differentiation occurs at shores, not at islands
# Irizarry RA, Ladd-Acosta C, Wen B, Wu Z, Montano C, Onyango P, Cui H, Gabo K, Rongione M, Webster M, Ji H, Potash JB, Sabunciyan S, Feinberg AP (2009). "The human colon cancer methylome shows similar hypo- and hypermethylation at conserved tissue-specific CpG island shores". Nature Genetics. 41 (2): 178–186

a <- table(feat$relation.to.island) / nrow(feat)
b <- table(feats$relation.to.island) / nrow(feats)
d1 <- as.data.frame(a)
d1$what <- "All"
d2 <- as.data.frame(b)
d2$what <- "Selected"
d <- rbind(d1, d2)

ggplot(d, aes(x=Var1, y=Freq)) +
geom_bar(stat="Identity", position="dodge", aes(fill=what))



# BMI
bmi <- scan("../data/bmi_cpg.txt", what="character")
table(bmi %in% cpglist)

table(table(top_sds$cpg[top_sds$cpg %in% bmi]))

cpgdat <- rbind(cpgdat,
	data.frame(cpg=bmi, source="BMI associated CpGs from Wilson et al 2017")
)


# Age

a <- read.csv("../data/AdditionalFile3.csv.gz")

age <- unique(as.character(a$CpGmarker[-1]))

table(age %in% cpglist)
table(table(top_sds$cpg[top_sds$cpg %in% age]))

cpgdat <- rbind(cpgdat,
	data.frame(cpg=age, source="Age associated CpGs used in Horvath predictor")
)

# Smoking (Illig)

load("../data/illig.RData")
smok <- Illig_data$cpgs

table(smok %in% cpglist)
table(table(top_sds$cpg[top_sds$cpg %in% smok]))

cpgdat <- rbind(cpgdat,
	data.frame(cpg=smok, source="Smoking associated CpGs from Zeilinger et al 2013")
)



## summarise and write

cpgsum <- dplyr::group_by(cpgdat, source) %>%
	dplyr::summarise(n=n())


save(cpgdat, file="../data/cpgdat.rdata")

write.csv(cpgsum, file="../data/cpgsum.csv")
