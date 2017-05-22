library(tidyverse)
library(meffil)
library(data.table)

# Obtain Zhou mapping file from here https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/45/4/10.1093_nar_gkw967/4/gkw967_Supplementary_Data.zip

# https://academic.oup.com/nar/article/45/4/e22/2290930/Comprehensive-characterization-annotation-and?rss=1#supplementary-data

# Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes

a <- fread("EPIC.anno.GRCh38.tsv", header=TRUE)
feat <- meffil.get.features()
b <- subset(a, probeID %in% feat$name & !MASK.general)
write.table(b$probeID, "../data/retain_from_zhou.txt", row=F, col=F, qu=F)

