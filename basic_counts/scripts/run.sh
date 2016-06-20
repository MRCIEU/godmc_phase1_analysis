
Rscript counts.R BSGS/results/05 BSGS/counts
Rscript counts.R LBC21/results/05 LBC21/counts
Rscript counts.R LBC36/results/05 LBC36/counts
Rscript counts.R ARIES/results/05 ARIES/counts
Rscript counts.R Leiden/results/05 Leiden/counts
Rscript counts.R SCZP1/results/05 SCZP1/counts

nom <- c("BSGS/counts.RData", "LBC21/counts.RData", "LBC36/counts.RData", "ARIES/counts.RData", "Leiden/counts.RData", "SCZP1/counts.RData")

l <- list()
for(i in 1:length(nom))
{
	load(nom[i])
	res$cohort <- strsplit(nom[i], split="/")[[1]][[1]]
	l[[i]] <- res
}

res <- do.call(rbind, l)

library(dplyr)

a <- group_by(subset(res, Var1 != "chr6"), pval, cohort) %>%
	dplyr::mutate(total = sum(count))

b <- group_by(a, pval, cohort, Var1) %>%
	dplyr::summarise(prop = count  / total)

library(meffil)

dat <- as.data.frame(table(meffil.featureset()$chromosome))

b <- merge(b, dat)

plot(prop ~ Freq, subset(b, pval==1e-8))

b$sex <- "Autosome"
b$sex[b$Var1 %in% c("chrX", "chrY")] <- "Sex"

b1 <- group_by(b, pval, cohort, sex) %>%
	dplyr::summarise(p = sum(prop))

plot(p ~ as.numeric(as.factor(sex)), b1)

ggplot(subset(b1, sex=="Sex"), aes(x=as.factor(pval), y=p)) +
geom_bar(stat="identity", aes(fill=cohort), position="dodge")

ggplot(subset(b, pval==1e-8), aes(x=Freq, y=prop)) +
geom_point(aes(colour=cohort, shape=sex), size=2) +
scale_colour_brewer(type="qual")



