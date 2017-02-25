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


neanderthal <- read.csv("../data/neanderthal-phewas-catalog.csv.gz")
neanderthal <- subset(neanderthal, !duplicated(snp), select=c(chr, bp, snp, maf_eur, maf_sas, maf_afr, maf_amr, maf_eas, genes))
nsn <- ucsc_get_position(neanderthal$snp)
neanderthal <- merge(nsn, neanderthal, by.y="snp", by.x="SNP")
save(neanderthal, file="../data/snps_neanderthal.rdata")
