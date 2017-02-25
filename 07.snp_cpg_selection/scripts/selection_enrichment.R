
### 

## Selection SNPs enrichment

test <- ucsc_get_position(unique(out$SNP))
test <- subset(test, !duplicated(name))
snp <- paste0(test$chrom, ":", test$chromEnd, ":SNP")
temp <- subset(dat, SNP %in% snp)
t.test(temp$fst[1:1000], dat$fst[sample(1:nrow(dat), 5000, replace=FALSE)])


data(aries_mqtl)
mqtl <- sample(unique(aries_mqtl$SNP), 5000, replace=FALSE)
mqtl2 <- ucsc_get_position(mqtl)
mqtl <- paste0(mqtl2$chrom, ":", mqtl2$chromEnd, ":SNP")
temp <- subset(dat, SNP %in% mqtl)
t.test(temp$fst[1:1000], dat$fst)

controls2 <- ucsc_get_position(controls)
controls <- paste0(controls2$chrom, ":", controls2$chromEnd, ":SNP")
temp <- subset(dat, SNP %in% controls)
t.test(temp$fst, dat$fst)


# It looks like mqtls, GWAS SNPs and hm3 genic control SNPs all have higher fst than average fst across the genome. 

# Just choose top fst SNPs that have maf > 0.1 in europeans

hist(selection$EUR)


