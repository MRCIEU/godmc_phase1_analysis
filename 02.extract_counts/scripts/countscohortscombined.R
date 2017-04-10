#path="/panfs/panasas01/shared-godmc/counts"
#l<-list.files(path,pattern="allcohorts.txt")

#g<-grep("cis",l)
#cis<-l[g]
#g<-grep("trans",l)
#trans<-l[g]
#r1.out<-data.frame()
#for (i in 1:length(cis)){
#r1<-read.table(cis[i],stringsAsFactors=F)
#spl<-strsplit(r1[,2],split="_")
#spl<-do.call("rbind",spl)
#r1<-data.frame(r1,spl)
#r1.out<-rbind(r1.out,r1)
#} 

#cohort summaries
#cd /panfs/panasas01/shared-godmc/cohort_summary

library("gdata")

path="/panfs/panasas01/shared-godmc/cohort_summary/"
coh<-read.table("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/CpG_SNP_FINAL/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts_noRAINE.txt",sep=" ",header=F)
coh<-coh[1:22,]
#mqtlcount.amcrae_BSGS.summary.txt

r.out<-data.frame()
for (i in 1:nrow(coh)){
r<-read.table(paste(path,"mqtlcount.",coh[i,1],"_",coh[i,2],".summary.txt",sep=""),header=T,stringsAsFactors=F)
r<-data.frame(coh[i,],r)
r.out<-rbind(r.out,r)
}

#extract N
path="/panfs/panasas01/shared-godmc/sftp/GoDMC"
#cohorts<-read.table("/panfs/panasas01/shared-godmc/scripts/cohorts.txt",sep=" ",header=F,stringsAsFactors=F)
cohorts<-coh
load(paste(path,cohorts[1,1],cohorts[1,2],"results/01/cohort_descriptives.RData",sep="/"))
df.out<-data.frame(names(cohort_summary))

for (i in 1:nrow(cohorts)){
#for (i in 1:2){
cat(i,"\n")
l<-list.files(paste(path,cohorts[i,1],cohorts[i,2],"results/01/",sep="/"),pattern="cohort_descriptives.RData")
if(length(l)==1){
load(paste(path,cohorts[i,1],cohorts[i,2],"results/01/cohort_descriptives.RData",sep="/"))
cohort_summary$predicted_cellcounts<-paste(cohort_summary$predicted_cellcounts,collapse=",")
cohort_summary$covariates<-paste(cohort_summary$covariates,collapse=",")
df <- data.frame(matrix(unlist(cohort_summary), nrow=length(cohort_summary), byrow=T),stringsAsFactors=FALSE)
df<-data.frame(names(cohort_summary),df)
df.out<-merge(df.out,df,by.x=1,by.y=1,all=T)
}
if(length(l)==0){
df<-data.frame(df.out[,1],NA)
df.out<-merge(df.out,df,by.x=1,by.y=1,all=T)

}
}
names(df.out)[-1]<-as.character(cohorts[,2])
#add TwinsUK
#df.out[6,"TwinsUK"]<-833

# Add in real N
real.N <- read.table("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/CpG_SNP_FINAL/repo/godmc_phase1_analysis/01.extract_sftp/data/cohort_samplesizes_noRAINE.txt")
studyN<-data.frame(names(df.out),N=t(df.out[6,]),N_snp=t(df.out[55,]),N_real=real.N[,3],cohort_real=real.N[,2])

#studyN<-data.frame(names(df.out),N=t(df.out[6,]),N_snp=t(df.out[55,]))

m<-match(r.out[,2],studyN[,1])
r.out<-data.frame(r.out,studyN[m,])
r.out <- drop.levels(r.out)

library(ggplot2)

#reorder_size <- function(x) {
#  factor(x, levels = names(sort(table(x))))
#}


#r.out$V2 <- as.factor(r.out$V2)
#r.out$X6<-as.factor(r.out$X6)
#a <- subset(r.out, !duplicated(V2), select=c(V2, X6))
#a <- a[order(a$X6),]
#levels(V2) <- a$V2

r.out2<-r.out
r.out<-r.out2
#r.out$V2<-paste(r.out$V2," (N=",r.out$X6,")",sep="")

df.o<-unique(data.frame(r.out$V2,r.out$N_real))
o<-order(as.numeric(as.character(df.o[,2])))
df.o<-df.o[o,]

r.out$V2 <- as.factor(r.out$V2)
r.out$N_real<-as.factor(r.out$N_real)
r.out$cisassoc_1000<-r.out$cisassoc/1000
r.out$cisassoc_1e6<-r.out$cisassoc/1000000
r.out$transassoc_1e6<-r.out$transassoc/1000000
r.out$V2<-factor(r.out$V2,levels=df.o[,1])
r.out$Cohort<-paste(r.out$V2," (N=",r.out$N_real,")",sep="")
r.out$Cohort<-gsub("Leiden_Longevity_Study","LLS",r.out$Cohort)

df.o<-unique(data.frame(r.out$Cohort,r.out$N_real))
o<-order(as.numeric(as.character(df.o[,2])))
df.o<-df.o[o,]

r.out$Cohort<-factor(r.out$Cohort,levels=df.o[,1])
m<-max(r.out$cisassoc_1e6)+1
p1<-ggplot(r.out,aes(x=reorder(pval,-pval),cisassoc_1e6,fill=Cohort)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of SNP-CpG pairs (cis)x1e6") +
scale_y_continuous(breaks=c(seq(0,m,2)), limits = c(0, m)) +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="cisassocbycohort.png",width=8,height=6)

m<-max(r.out$transassoc_1e6)+1
p1<-ggplot(r.out,aes(x=reorder(pval,-pval),transassoc_1e6,fill=Cohort)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of SNP-CpG pairs (trans)x1e6") +
scale_y_continuous(breaks=c(seq(0,m,2)), limits = c(0, m)) +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="transassocbycohort.png",width=8,height=6)

r.out$cisCpGs_1000<-r.out$cisCpGs/1000
r.out$transCpGs_1000<-r.out$transCpGs/1000

p1<-ggplot(r.out,aes(x=reorder(pval,-pval),cisCpGs_1000,fill=Cohort)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of CpGs (cis)x1000") +
scale_y_continuous(breaks=c(seq(0,max(r.out$cisCpGs_1000)+1,10)), limits = c(0, max(r.out$cisCpGs_1000)+1)) +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="cisCpGsbycohort.png",width=8,height=6)

p1<-ggplot(r.out,aes(x=reorder(pval,-pval),transCpGs_1000,fill=Cohort)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of CpGs (trans)x1000") +
scale_y_continuous(breaks=c(seq(0,max(r.out$transCpGs_1000)+1,50)), limits = c(0, max(r.out$transCpGs_1000)+1)) +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="transCpGsbycohort.png",width=8,height=6)


r.out_n<-unique(r.out[,c(1,2,10,13,12,18)])
r.out_n$N_real<-as.numeric(as.character(r.out_n$N_real))
r.out_n$X55<-as.numeric(as.character(r.out_n$X55))

p1<-ggplot(r.out_n,aes(x=N_real,y=X55))+
geom_point(aes(colour=factor(r.out_n$Cohort),size=r.out_n$N_real)) +
xlab("Cohort N") +
ylab("N SNPs") +
#scale_y_continuous(breaks=c(seq(0,max(r.out$transCpGs_1000)+1,50)), limits = c(0, max(r.out$transCpGs_1000)+1)) +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="SNPsbyNcohort.png",width=8,height=6)






