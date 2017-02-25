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


hm3_genic <- scan("../data/genicsnps.txt.gz", what="character")
hm3_genic <- sample(hm3_genic, 2000, replace=FALSE)
table(hm3_genic %in% gcc$SNP)
hm3_genic <- ucsc_get_position(hm3_genic)
hm3_genic <- hm3_genic[! hm3_genic$SNP %in% gcc$SNP, ][1:1000, ]

write.table(hm3_genic, file = "../data/snps_control.txt", row=F, col=F, qu=F)
