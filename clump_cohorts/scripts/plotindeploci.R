library(ggplot2)
library(meffil)

y<-meffil.get.features("450k")
ncpg<-length(unique(y$name))


####
n_independent_regions <- 1000000

#3 billion basepairs residing in 23 pairs of chromosomes
n_bases <- 3000000000

#SNP-CpG distance is 1 Mb 
cis_window <- 1000000

n_independent_regions_cis <- n_independent_regions / (n_bases / cis_window)
n_independent_regions_trans <- n_independent_regions - n_independent_regions_cis

#number of analysed CpGs - all probe on 450k array
ncpg <- ncpg

ntest_cis <- ncpg * n_independent_regions_cis
ntest_trans <- ncpg * n_independent_regions_trans

# we used pval threshold of 1e-05
exp_cis <- ntest_cis * 1e-5
exp_trans <- ntest_trans * 1e-5
exp_cis
exp_trans

#

path="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/clump_cohorts/data"
l<-list.files(path=path,pattern=".numberofindependentloci.Robj")

cohorts<-read.table("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts.txt",sep=" ",header=F)
cohorts<-cohorts[1:length(l),]


indexsnps<-data.frame()
for (i in 1:length(l)){

load(paste(path,"/",cohorts[i,1],"_",cohorts[i,2],".numberofindependentloci.Robj",sep=""))

noindexsnps<-sum(r.out$snpno,na.rm=T)
nocisindexsnps<-sum(r.out$cis,na.rm=T)
notransindexsnps<-sum(r.out$trans,na.rm=T)
df<-data.frame(user=cohorts[i,1],cohort=cohorts[i,2],noindexsnps,nocisindexsnps,notransindexsnps)
indexsnps<-rbind(indexsnps,df)
}

#extract N and number of SNPs
path="/panfs/panasas01/shared-godmc/sftp/GoDMC"
#cohorts<-read.table("/panfs/panasas01/shared-godmc/scripts/cohorts.txt",sep=" ",header=F,stringsAsFactors=F)

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
studyN<-data.frame(names(df.out),N=t(df.out[6,]),N_snp=t(df.out[55,]))

m<-match(indexsnps$cohort,studyN[,1])
r.out<-data.frame(indexsnps,studyN[m,])



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
r.out$cohort<-gsub("Leiden_Longevity_Study","LLS",r.out$cohort)

df.o<-unique(data.frame(r.out$cohort,r.out$X6))
o<-order(as.numeric(as.character(r.out$X6)))
r.out<-r.out[o,]

r.out$cohort <- as.factor(r.out$cohort)
r.out$X6<-as.factor(r.out$X6)
r.out$nocisindexsnps_1000<-r.out$nocisindexsnps/1000
r.out$nocisindexsnps_1e6<-r.out$nocisindexsnps/1000000
r.out$notransindexsnps_1e6<-r.out$notransindexsnps/1000000
r.out$noindexsnps_1e6<-r.out$noindexsnps/1000000

r.out$cohort<-factor(r.out$cohort,levels=df.o[,1])
r.out$Cohort<-paste(r.out$cohort," (N=",r.out$X6,")",sep="")
r.out$Cohort<-gsub("Leiden_Longevity_Study","LLS",r.out$Cohort)


df.o<-unique(data.frame(r.out$Cohort,r.out$X6))
o<-order(as.numeric(as.character(df.o[,2])))
df.o<-df.o[o,]

r.out$Cohort<-factor(r.out$Cohort,levels=df.o[,1])
m<-max(r.out$nocisindexsnps_1000)+1
p1<-ggplot(r.out,aes(x=cohort,nocisindexsnps_1000,fill=Cohort)) +
geom_bar(stat="identity",position="dodge") +
xlab("Cohort") +
ylab("number of indexSNPs (cis) x1e3") +
geom_hline(yintercept = exp_cis/1000) +
scale_y_continuous(breaks=c(seq(0,m,10)), limits = c(0, m)) +
theme(axis.text.x = element_text(face = "bold",angle=90,hjust=1))
ggsave(p1,file="cisindexSNPbycohort.png",width=8,height=6)

r.out$Cohort<-factor(r.out$Cohort,levels=df.o[,1])
m<-max(r.out$notransindexsnps_1e6)+1
p1<-ggplot(r.out,aes(x=cohort,notransindexsnps_1e6,fill=Cohort)) +
geom_bar(stat="identity",position="dodge") +
xlab("Cohort") +
ylab("number of indexSNPs (trans) x1e6") +
geom_hline(yintercept = exp_trans/1000000) +
scale_y_continuous(breaks=c(seq(0,m,10)), limits = c(0, m)) +
theme(axis.text.x = element_text(face = "bold",angle=90,hjust=1))
ggsave(p1,file="transindexSNPbycohort.png",width=8,height=6)

r.out$Cohort<-factor(r.out$Cohort,levels=df.o[,1])
m<-max(r.out$noindexsnps_1e6)+1
p1<-ggplot(r.out,aes(x=cohort,noindexsnps_1e6,fill=Cohort)) +
geom_bar(stat="identity",position="dodge") +
xlab("Cohort") +
ylab("number of indexSNPs x1e6") +
scale_y_continuous(breaks=c(seq(0,m,10)), limits = c(0, m)) +
theme(axis.text.x = element_text(face = "bold",angle=90,hjust=1))
ggsave(p1,file="indexSNPbycohort.png",width=8,height=6)

