library(meffil)
args <- (commandArgs(TRUE));
user <- as.character(args[1]);
cohort <- as.character(args[2]);


y<-meffil.featureset("450k")

users<-user
cohorts=cohort

#users=c("epzjlm","ehannon","rluijk","ecarnero")
#cohorts=c("ARIES","Phase1_SCZ","Leiden_Longevity_Study","TwinsUK")

for (k in 1:length(users)){
#k=1
path=paste("/panfs/panasas01/shared-godmc/sftp/GoDMC/",users[k],"/",cohorts[k],"/results/05/",sep="")

l<-list.files(path=path,pattern="RData")

eqtl<-data.frame()

pvals<-c(1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13)
#pvals<-c(1e-13)
eqtl.count<-matrix(nrow=length(pvals)*length(l),ncol=7)
s<-seq(0,length(l)*length(pvals),length(pvals))

probe.trans.all<- probe.cis.all <- lapply(1:length(pvals), function(i) character(0))
assoc.trans.all<- assoc.cis.all <- lapply(1:length(pvals), function(i) character(0))

for (i in (1:length(l))){
#for (i in (384:500)){

load(paste(path,"res.",i,".RData",sep=""))
eqtl<-me$all$eqtls

#names(eqtl)
#[1] "snps"      "gene"      "statistic" "pvalue"    "FDR"       "beta"     

#ID	MARKER	RSID	CHR	POS	CPG	CPGchr	CPGpos	BETA	tstat	P	FDR	Trans	EA	NEA	EAF	SE
#SNP	CpG	MARKER	BETA	SE	P	CHR	SNP	EA	NEA	EAF	N

assoc<-paste(eqtl$snps,eqtl$gene,sep="_")
m<-match(me$all$eqtls$gene,y$name)
spl<-strsplit(as.character(eqtl$snps),split=":")
spl<-do.call("rbind",spl)

#    [,1]   [,2]      [,3]   
#[1,] "chr1" "2578220" "SNP"  
#[2,] "chr1" "2694496" "SNP"  


eqtl<-data.frame(assoc,spl,me$all$eqtls,y[m,c("name","chromosome","position","relation.to.island")])
eqtl$X1<-gsub("chr","",eqtl$X1)
eqtl$chromosome<-gsub("chr","",eqtl$chromosome)
eqtl$chromosome<-gsub("X","23",eqtl$chromosome)
eqtl$chromosome<-gsub("Y","24",eqtl$chromosome)
eqtl$X2<-as.numeric(as.character(eqtl$X2))
eqtl$SE = eqtl$beta / eqtl$statistic
w<-which(names(eqtl)%in%c("X3","statistic","FDR","name"))
eqtl<-eqtl[,-w]
eqtl2<-data.frame(SNP=eqtl$snps,CpG=eqtl$gene,MARKER=eqtl$assoc,BETA=eqtl$beta,SE=eqtl$SE,P=eqtl$pvalue)

chr<-which(eqtl$X1==eqtl$chromosome)
eqtl$cischr<-ifelse(eqtl$X1==eqtl$chromosome,1,0)
eqtl$cis<-ifelse(eqtl$cischr==1 & abs(eqtl$X2-eqtl$position)<1000000,"cis","trans")
#trans chr
#eqtl$transchr<-ifelse(eqtl$cischr==0,"trans","cis")

write.table(eqtl2,paste("/panfs/panasas01/shared-godmc/meta-analysis/inputfiles/",users[k],"_",cohorts[k],".",i,".gwama.formatted.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)

probe.cis<-list()
probe.trans<-list()
assoc.cis<-list()
assoc.trans<-list()

for (p in 1:length(pvals)){
s1<-s[i]

cat(pvals[p],"\n")
cat(s1+p,"\n")
eqtl<-eqtl[which(eqtl$pvalue<as.numeric(pvals[p])),]

eqtl.count[(s1+p),1]<-pvals[p]
eqtl.count[(s1+p),2]<-as.character(table(eqtl$cis)[1])
eqtl.count[(s1+p),3]<-as.character(table(eqtl$cis)[2])
eqtl.count[(s1+p),4]<-length(unique(eqtl[eqtl$cis=="cis","snps"]))
eqtl.count[(s1+p),5]<-length(unique(eqtl[eqtl$cis=="trans","snps"]))
eqtl.count[(s1+p),6]<-length(unique(eqtl[eqtl$cis=="cis","gene"]))
eqtl.count[(s1+p),7]<-length(unique(eqtl[eqtl$cis=="trans","gene"]))

probe.cis[[p]]<-unique(as.character(eqtl[eqtl$cis=="cis","gene"]))
probe.trans[[p]]<-unique(as.character(eqtl[eqtl$cis=="trans","gene"]))

assoc.cis[[p]]<-unique(as.character(eqtl[eqtl$cis=="cis","assoc"]))
assoc.trans[[p]]<-unique(as.character(eqtl[eqtl$cis=="trans","assoc"]))

write.table(assoc.cis[[p]],paste("/panfs/panasas01/shared-godmc/counts/",users[k],"_",cohorts[k],"/cis.assoc.",i,".",users[k],"_",cohorts[k],".",pvals[p],".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
write.table(assoc.trans[[p]],paste("/panfs/panasas01/shared-godmc/counts/",users[k],"_",cohorts[k],"/trans.assoc.",i,".",users[k],"_",cohorts[k],".",pvals[p],".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

}


#probe.cis.all<-lapply(1:length(probe.cis.all), function(i) c(probe.cis.all[[i]], probe.cis[[i]]))
probe.cis.all<-mapply(c,probe.cis.all,probe.cis)
probe.trans.all<-mapply(c,probe.trans.all,probe.trans)

assoc.cis.all<-mapply(c,assoc.cis.all,assoc.cis)
assoc.trans.all<-mapply(c,assoc.trans.all,assoc.trans)

probe.cis.all <- lapply(probe.cis.all, unique)
cis.probes<-unlist(sapply(probe.cis.all, length))

probe.trans.all <- lapply(probe.trans.all, unique)
trans.probes<-unlist(sapply(probe.trans.all, length))

assoc.cis.all <- lapply(assoc.cis.all, unique)
assoc.trans.all <- lapply(assoc.trans.all, unique)

}
write.table(eqtl.count,paste("/panfs/panasas01/shared-godmc/cohort_summary/mqtlcount.",users[k],"_",cohorts[k],".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

#save(assoc.cis.all,file=paste("/panfs/panasas01/shared-godmc/counts/cis.assoc.",users[k],".Robj",sep=""))
#save(assoc.trans.all,file=paste("/panfs/panasas01/shared-godmc/counts/trans.assoc.",users[k],".Robj",sep=""))

###

eqtl.count<-read.table(paste("/panfs/panasas01/shared-godmc/cohort_summary/mqtlcount.",users[k],"_",cohorts[k],".txt",sep=""),header=T)
eqtl.count<-data.frame(eqtl.count)
for (i in 1:7){eqtl.count[,i]<-as.numeric(eqtl.count[,i])}

out<-data.frame()
for (p in 1:length(pvals)){ 
eqtl.count.p<-eqtl.count[eqtl.count[,1]==pvals[p],-1]
a<-apply(eqtl.count.p,2,function(x) sum(x,na.rm=T))
out<-rbind(out,a)}

out<-data.frame(pvals,out)
names(out)<-c("pval","cisassoc","transassoc","cisSNPs","transSNPs","cisCpGs","transCpGs")
out$cisCpGs<-cis.probes
out$transCpGs<-trans.probes

write.table(out,paste("/panfs/panasas01/shared-godmc/cohort_summary/mqtlcount.",users[k],"_",cohorts[k],".summary.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
}


