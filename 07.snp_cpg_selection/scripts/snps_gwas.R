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


ao <- available_outcomes()
a <- subset(ao,
	priority == 1 &
	category %in% c("Disease", "Risk factor") &
	population %in% c("European", "Mixed") & 
	! id %in% c(981, 1026, 1080
		# , 1089, 1090, 1031, 1085, 27, 1083, 276, 278, 1059, 1060, 991, 1074, 9, 10, 12, 980, 1013, 294, 982, 965, 966, 985, 986, 280, 286, 1025, 118, 811, 281, 283, 832, 1076, 1077, 1078, 1071, 1072, 967, 815, 32, 968, 755, 837
	) &
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
mrbase$trait <- strsplit(mrbase$exposure, split="\\|") %>% sapply(function(x) x[1]) %>% gsub(" $", "", .)


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

save(mrbase, file="")

data(gwas_catalog)
g <- filter(gwas_catalog, pval < 5e-8) %>%
	group_by(Phenotype_simple) %>%
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

gwascat <- format_data(g)
gwascat$data_source.exposure <- "GWAS catalog"

gwas <- bind_rows(gwascat, out)
snpid <- ucsc_get_position(unique(gwas$SNP))
gwas <- merge(gwas, snpid)

save(gwas, file="../data/snps_gwas.rdata")

