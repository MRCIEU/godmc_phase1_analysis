p<-read.table("SOMA_pQTL.txt",he=T)
p$POS<-gsub(",","",p$POS)
p$POS<-as.numeric(as.character(p$POS))

n1<-nchar(as.character(p$EA))
n2<-nchar(as.character(p$NEA))
p<-data.frame(p,n1,n2)
l<-which(p$n1=="1"&p$n2=="1")
p$id[l]<-paste("chr",p$CHR[l],":",p$POS[l],":","SNP",sep="")
p$id[-l]<-paste("chr",p$CHR[-l],":",p$POS[-l],":","INDEL",sep="")
save(p,file="SOMA_pQTL.RData")