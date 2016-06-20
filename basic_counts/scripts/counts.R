arguments <- commandArgs(T)
indr <- arguments[1]
otnm <- arguments[2]
fn1 <- paste0(otnm, ".RData")
fn2 <- paste0(otnm, ".pdf")


library(meffil)
library(dplyr)
library(ggplot2)

files <- dir(indr) %>% grep(".RData", ., value=TRUE) %>% paste0(indr, "/", .)
featureset <- subset(meffil.featureset(), select=c(name, chromosome))
pvals <- 10^-(5:14)

count_res <- function(filename, indr, featureset, pvals)
{
	message(filename)
	load(filename)
	res1 <- merge(me$all$eqtl, featureset, by.x="gene", by.y="name")

	a <- lapply(pvals, function(x) 
		{
			d <- subset(res1, pvalue < x)
			d <- data.frame(table(d$chromosome))
			d$pval <- x
			return(d)
		}) %>% bind_rows

	return(a)
}

res <- mclapply(files, function(x) count_res(x, indr, featureset, pvals), mc.cores=16)

res <- bind_rows(res) %>% group_by(Var1, pval) %>% summarise(count=sum(Freq))


p1 <- ggplot(res, aes(x=Var1, y=count)) + geom_bar(stat="identity", aes(fill=as.factor(pval)), position="dodge") + labs(x="Chromsome")

save(res, file=fn1)
ggsave(p1, filename=fn2, height=6, width=12)

