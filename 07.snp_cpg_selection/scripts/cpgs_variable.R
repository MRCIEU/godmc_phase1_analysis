library(tidyverse)
library(meffil)
library(data.table)

feat <- meffil.get.features()
#486425
xy<-which(feat$chromosome%in%c("chrX","chrY"))
probes_xy<-as.character(feat[xy,"name"])
#11648
feat_noxy<-feat[-xy,]
#474777

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
	w<-which(meth_summary$cpg%in%probes_xy)
	if(length(w)>0){meth_summary<-meth_summary[-w,]}
	meth_summary$cohort <- info$cohort[i]

	l[[i]] <- meth_summary
}

dat <- bind_rows(l)

means <- spread(subset(dat, select=c(cpg, cohort, mean)), key=cohort, value=mean)
sds <- spread(subset(dat, select=c(cpg, cohort, sd)), key=cohort, value=sd)

cormeans <- cor(means[,-c(1)], use="pair")
corsds <- cor(sds[,-c(1)], use="pair")


# load("../data/filtered_probe_list.RData")
# retain <- scan("../data/retain_from_Naeem_multimap.txt", what="character")
# mafs <- read.table("../data/SNPatCpG_AFs.txt.gz", he=T, stringsAsFactors=FALSE)
# retaincpg <- unique(c(probeinfo$TargetID, retain))

# mafs2 <- strsplit(mafs$EUR_AF, split=",")
# rem <- sapply(mafs2, function(x) any(x > 0.01 & x < 0.99))
# retaincpg <- retaincpg[! retaincpg %in% mafs$name[rem]]


retaincpg <- scan("../data/retain_from_zhou.txt", what="character")
#435391

#exclusion probes from TwinsUK
excl<-read.table("../data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509

dat_clean <- subset(dat, cpg %in% retaincpg)
dim(dat_clean)
top_sds <- group_by(dat_clean, cohort) %>%
	do({
		x <- .
		thresh <- quantile(x$sd, 0.1)
		return(subset(x, sd >= thresh, select=c(cpg, sd)))
	})

tab_cpg <- table(top_sds$cpg)
table(tab_cpg)
cpglist <- names(tab_cpg)[tab_cpg >= 20]
length(cpglist)

##

cpgdat <- data.frame(cpg = cpglist, source = "320k most variable sites amongst all 23 GoDMC phase 1 cohorts", stringsAsFactors=FALSE)

##




# Randomly sampled from each chromosome

feats <- subset(feat_noxy, name %in% cpglist)

a <- table(feat_noxy$chromosome) / nrow(feat_noxy)
b <- table(feats$chromosome) / nrow(feats)

pdf("../data/cpgsselection.pdf",height=6,width=6)
plot(as.numeric(a), as.numeric(b), xlab='Proportion of CpGs on chromosome', ylab="Proportion of selected CpGs on chromosome")


# Selected has more from open sea and shores, less from islands

# Most tissue / cancer differentiation occurs at shores, not at islands
# Irizarry RA, Ladd-Acosta C, Wen B, Wu Z, Montano C, Onyango P, Cui H, Gabo K, Rongione M, Webster M, Ji H, Potash JB, Sabunciyan S, Feinberg AP (2009). "The human colon cancer methylome shows similar hypo- and hypermethylation at conserved tissue-specific CpG island shores". Nature Genetics. 41 (2): 178–186

a <- table(feat_noxy$relation.to.island) / nrow(feat_noxy)
b <- table(feats$relation.to.island) / nrow(feats)
d1 <- as.data.frame(a)
d1$what <- "All"
d2 <- as.data.frame(b)
d2$what <- "Selected"
d <- rbind(d1, d2)

ggplot(d, aes(x=Var1, y=Freq)) +
geom_bar(stat="Identity", position="dodge", aes(fill=what))

dev.off()

# heritable probes
# filtered on h2_total>0.5
load("../data/vanDongen_h2gt50_probes.RData")
h2probes<-as.character(r$cgid)

table(h2probes %in% cpglist)
table(h2probes %in% retaincpg)

table(table(top_sds$cpg[top_sds$cpg %in% h2probes]))

cpgdat <- rbind(cpgdat,
	data.frame(cpg=h2probes, source="heritable CpGs (h2>0.5) from van Dongen et al 2016")
)

#For this we used 364 twin pairs adjusting for age, BMI, cell heterogeneity, smoking, alcohol, and technical covs (plate and position). 
h<-read.table("../data/TwinsUK_330MZ_34DZ_h2.txt",sep="\t",he=F)
h<-h[which(h$V2>0.5),]
#28618
h<-as.character(h$V1)

table(h %in% cpglist)
table(h %in% retaincpg)

table(table(top_sds$cpg[top_sds$cpg %in% h]))

cpgdat <- rbind(cpgdat,
	data.frame(cpg=h, source="heritable CpGs (h2>0.5) from TwinsUK")
)


# BMI
bmi <- scan("../data/bmi_cpg.txt", what="character")
table(bmi %in% cpglist)
table(bmi %in% retaincpg)

table(table(top_sds$cpg[top_sds$cpg %in% bmi]))

cpgdat <- rbind(cpgdat,
	data.frame(cpg=bmi, source="BMI associated CpGs from Wahl et al 2017")
)


# Age

a <- read.csv("../data/AdditionalFile3.csv.gz")

age <- unique(as.character(a$CpGmarker[-1]))

table(age %in% cpglist)
table(age %in% retaincpg)
table(table(top_sds$cpg[top_sds$cpg %in% age]))

cpgdat <- rbind(cpgdat,
	data.frame(cpg=age, source="Age associated CpGs used in Horvath predictor")
)

# Smoking (Illig)

load("../data/illig.RData")
smok <- Illig_data$cpgs

load("../data/joehanes.rdata")
smok <- joehanes$Probe.ID

table(smok %in% cpglist)
table(table(top_sds$cpg[top_sds$cpg %in% smok]))
table(smok %in% retaincpg)
cpgdat <- rbind(cpgdat,
	data.frame(cpg=smok, source="Smoking associated CpGs from Joehanes et al")
)


# scz

a1 <- read.table("../data/scz_cpg1.txt", sep="\t", stringsAsFactors=FALSE)
a2 <- read.table("../data/scz_cpg2.txt", sep=" ", stringsAsFactors=FALSE)

table(a1$V1 %in% cpglist)
table(a1$V1 %in% retaincpg)
table(a2$V1 %in% cpglist)
table(a2$V1 %in% retaincpg)

cpgdat <- rbind(cpgdat,
	data.frame(cpg=a1$V1, source="Schizophrenia, Hannon et al 2016"),
	data.frame(cpg=a2$V1, source="Schizophrenia, Montana et al 2016")
)

## summarise and write

cpgsum <- dplyr::group_by(cpgdat, source) %>%
	dplyr::summarise(n=n())
cpgsum <- rbind(cpgsum, data.frame(source="Total unique CpGs", n=length(unique(cpgdat$cpg))))

save(cpgdat, file="../data/cpgs_variable.rdata")

write.csv(cpgsum, file="../data/cpgsum.csv")
