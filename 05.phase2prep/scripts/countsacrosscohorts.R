###make plots for counts across cohorts

#/panfs/panasas01/shared-godmc/counts_summary
library(ggplot2)
path="/panfs/panasas01/shared-godmc/counts_2017/combined"
cohorts<-read.table("/panfs/panasas01/shared-godmc/scripts/cohorts_noRAINE.txt",sep=" ",header=F)
cohorts<-cohorts[1:22,]

counts<-data.frame()
for (i in 1:nrow(cohorts)){
r<-read.table(paste(path,"/counts.allcohorts.combined.",i,".txt",sep=""),header=F,stringsAsFactors=F) #rows are 9 pval threshold*7
cat(nrow(r),"\n")

w<-which(r[,1]%in%"Pvalue")
w2<-nrow(r)/length(w)
s<-w+(w2-1)

r.out<-data.frame()
for (j in 1:length(w)){
r2<-r[w[j]:s[j],]
lab<-r2[,1]
r2<-t(r2[,2])
colnames(r2)<-lab
r.out<-rbind(r.out,r2)
}
r.out<-data.frame(nocohorts=i,r.out)
counts<-rbind(counts,r.out)
}

counts$Overlapping_cis_SNP.CpG_pairs.N_1e6<-counts$Overlapping_cis_SNP.CpG_pairs.N./1e6
counts$Overlapping_trans_SNP.CpG_pairs.N_1e6<-counts$Overlapping_trans_SNP.CpG_pairs.N./1e6
counts$nocohorts<-as.factor(counts$nocohorts)

p1<-ggplot(counts,aes(x=reorder(Pvalue,-Pvalue),Overlapping_cis_SNP.CpG_pairs.N_1e6,fill=nocohorts)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of overlapping SNP-CpG pairs (cis)x1e6") +
scale_y_continuous(breaks=c(seq(0,max(counts$Overlapping_cis_SNP.CpG_pairs.N_1e6)+1,5)), limits = c(0, max(counts$Overlapping_cis_SNP.CpG_pairs.N_1e6)+1)) +
scale_fill_discrete(name="No of cohorts") +
theme(axis.text.x = element_text(face = "bold"))

ggsave(p1,file="/panfs/panasas01/shared-godmc/counts_summary_2017/cisassocacrosscohort.png",height=6,width=8)

p1<-ggplot(counts,aes(x=reorder(Pvalue,-Pvalue),Overlapping_trans_SNP.CpG_pairs.N_1e6,fill=nocohorts)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of overlapping SNP-CpG pairs (trans)x1e6") +
scale_y_continuous(breaks=c(seq(0,max(counts$Overlapping_trans_SNP.CpG_pairs.N_1e6)+1,50)), limits = c(0, max(counts$Overlapping_trans_SNP.CpG_pairs.N_1e6)+1)) +
scale_fill_discrete(name="No of cohorts") +
theme(axis.text.x = element_text(face = "bold"))

ggsave(p1,file="/panfs/panasas01/shared-godmc/counts_summary_2017/transassocacrosscohort.png",height=6,width=8)

p1<-ggplot(counts,aes(x=reorder(Pvalue,-Pvalue),Overlapping_cis_SNP.CpG_pairs...,fill=nocohorts)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of overlapping SNP-CpG pairs (%) (cis) ") +
scale_y_continuous(breaks=c(seq(0,100,10)), limits = c(0, 100)) +
scale_fill_discrete(name="No of cohorts") +
theme(axis.text.x = element_text(face = "bold"))

ggsave(p1,file="/panfs/panasas01/shared-godmc/counts_summary_2017/cisassocacrosscohortsperc.png",height=6,width=8)

p1<-ggplot(counts,aes(x=reorder(Pvalue,-Pvalue),Overlapping_trans_SNP.CpG_pairs...,fill=nocohorts)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of overlapping SNP-CpG pairs (%) (trans) ") +
scale_y_continuous(breaks=c(seq(0,100,10)), limits = c(0,100)) +
scale_fill_discrete(name="No of cohorts") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="/panfs/panasas01/shared-godmc/counts_summary_2017/transassocacrosscohortsperc.png",height=6,width=8)

#counts cumulative
counts.c<-data.frame()
pvals<-unique(counts$Pvalue)
for (i in 1:length(pvals)){
counts.p<-counts[which(counts$Pvalue==pvals[i]),]
o<-order(counts.p$nocohorts,decreasing=T)
counts.p<-counts.p[o,]
c1<-cumsum(counts.p[,3])
c2<-cumsum(counts.p[,4])
df<-data.frame(nocohorts=nrow(counts.p):1,Pvalue=counts.p$Pvalue,Overlapping_cis_SNP.CpG_pairs.N=c1,Overlapping_trans_SNP.CpG_pairs.N=c2,counts.p[,5:6])
o<-order(df$nocohorts)
df<-df[o,]

counts.c<-rbind(counts.c,df)
}

counts.c$Overlapping_cis_SNP.CpG_pairs.N_1e6<-counts.c$Overlapping_cis_SNP.CpG_pairs.N/1e6
counts.c$Overlapping_trans_SNP.CpG_pairs.N_1e6<-counts.c$Overlapping_trans_SNP.CpG_pairs.N/1e6
counts.c$nocohorts<-as.factor(counts.c$nocohorts)
counts.c$nocohorts.lab<-paste(">",counts.c$nocohorts)

p1<-ggplot(counts.c,aes(x=reorder(Pvalue,-Pvalue),Overlapping_cis_SNP.CpG_pairs.N_1e6,fill=nocohorts)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of overlapping SNP-CpG pairs (cis)x1e6") +
scale_y_continuous(breaks=c(seq(0,max(counts.c$Overlapping_cis_SNP.CpG_pairs.N_1e6)+1,2)), limits = c(0, max(counts.c$Overlapping_cis_SNP.CpG_pairs.N_1e6)+1)) +
scale_fill_discrete(name="No of cohorts") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="/panfs/panasas01/shared-godmc/counts_summary_2017/cisassocacrosscohorts_cumulative.png",height=6,width=8)

p1<-ggplot(counts.c,aes(x=reorder(Pvalue,-Pvalue),Overlapping_trans_SNP.CpG_pairs.N_1e6,fill=nocohorts)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of overlapping SNP-CpG pairs (trans)x1e6") +
scale_y_continuous(breaks=c(seq(0,max(counts.c$Overlapping_trans_SNP.CpG_pairs.N_1e6)+1,8)), limits = c(0, max(counts.c$Overlapping_trans_SNP.CpG_pairs.N_1e6)+1)) +
scale_fill_discrete(name="No of cohorts") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="/panfs/panasas01/shared-godmc/counts_summary_2017/transassocacrosscohorts_cumulative.png",height=6,width=8)

counts.c2<-counts.c[counts.c$nocohorts!=1,]

p1<-ggplot(counts.c2,aes(x=reorder(Pvalue,-Pvalue),Overlapping_cis_SNP.CpG_pairs.N_1e6,fill=nocohorts)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of overlapping SNP-CpG pairs (cis)x1e6") +
scale_y_continuous(breaks=c(seq(0,max(counts.c2$Overlapping_cis_SNP.CpG_pairs.N_1e6)+1,2)), limits = c(0, max(counts.c2$Overlapping_cis_SNP.CpG_pairs.N_1e6)+1)) +
scale_fill_discrete(name="No of cohorts") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="/panfs/panasas01/shared-godmc/counts_summary_2017/cisassocacrosscohorts_cumulative_atleast2.png",height=6,width=8)

p1<-ggplot(counts.c2,aes(x=reorder(Pvalue,-Pvalue),Overlapping_trans_SNP.CpG_pairs.N_1e6,fill=nocohorts)) +
geom_bar(stat="identity",position="dodge") +
xlab("pvalue threshold") +
ylab("number of overlapping SNP-CpG pairs (trans)x1e6") +
scale_y_continuous(breaks=c(seq(0,max(counts.c2$Overlapping_trans_SNP.CpG_pairs.N_1e6)+1,2)), limits = c(0, max(counts.c2$Overlapping_trans_SNP.CpG_pairs.N_1e6)+1)) +
scale_fill_discrete(name="No of cohorts") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="/panfs/panasas01/shared-godmc/counts_summary_2017/transassocacrosscohorts_cumulative_atleast2.png",height=6,width=8)


#any association
#p<-length(unique(counts$Pvalue))
#df<-data.frame(nocohorts=rep("Any",p),Pvalue=counts$Pvalue[1:p],Overlapping_cis_SNP.CpG_pairs.N.=counts$All_cis_SNP.CpG_pairs.N.[1:p],Overlapping_trans_SNP.CpG_pairs.N.=counts$All_trans_SNP.CpG_pairs.N.[1:p],All_cis_SNP.CpG_pairs.N.=counts$All_cis_SNP.CpG_pairs.N.[1:p],All_trans_SNP.CpG_pairs.N.=counts$All_trans_SNP.CpG_pairs.N.[1:p],Overlapping_cis_SNP.CpG_pairs...=rep(100,p),Overlapping_trans_SNP.CpG_pairs...=rep(100,p))
#counts$nocohorts<-as.factor(counts$nocohorts)
#counts<-rbind(df,counts)





