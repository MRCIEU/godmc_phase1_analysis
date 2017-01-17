library(meffil)
library(data.table)

cis_distance <- 1000000
cpg_info <- meffil.featureset("450k")
file_list <- file.path("~/sandpit/splitfiles", list.files("~/sandpit/splitfiles/", pattern="gz"))

for(i in 1:length(file_list))
{
	message(i, " of ", length(file_list))
	l <- fread(paste0("zcat ", file_list[i]))
	n <- nrow(l)
	l$index <- 1:nrow(l)
	l <- merge(l, subset(cpg_info, select=c(name, chromosome, position)), by.x="V3", by.y="name")
	l <- l[order(l$index), ]
	stopifnot(n == nrow(l))
	temp <- do.call(rbind, strsplit(l$V2, split=":"))
	l$chr <- temp[,1]
	l$pos <- as.numeric(temp[,2])
	l$cis <- l$chr == l$chromosome & (abs(l$pos - l$position) < cis_distance)
	l$cis <- as.character(l$cis)
	l$cis[l$cis == "TRUE"] <- "c"
	l$cis[l$cis == "FALSE"] <- "t"
	l <- subset(l, select=c(V1, V2, V3, cis))
	outfile <- gzfile(file_list[i], "w")
	write.table(l, outfile, row=FALSE, col=FALSE, qu=FALSE)
	close(outfile)
}
