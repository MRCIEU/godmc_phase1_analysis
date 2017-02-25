library(TwoSampleMR)
library(MRInstruments)
library(dplyr)



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

expression <- format_gtex_eqtl(subset(gtex_eqtl, tissue == "Whole Blood" & pval < 5e-10))
snpid <- ucsc_get_position(unique(expression$SNP))
expression <- merge(expression, snpid)

save(expression, file="../data/snps_expression.rdata")
