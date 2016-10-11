library(data.table)
path="/panfs/panasas01/shared-godmc/meta-analysis/betacomparison/"
l<-list.files(path=path,pattern=".gwama.formatted.chr20.txt")
l2<-gsub(".gwama.formatted.chr20.txt","",l)
l2<-gsub("Leiden_Longevity_Study","LLS",l2)
l2<-gsub("Phase1_SCZ","Phase1SCZ",l2)
l2<-gsub("Phase2_SCZ","Phase2SCZ",l2)

spl<-strsplit(l2,split="_")
spl<-do.call("rbind",spl)

ref<-read.table("~/repo/godmc/resources/genetics/1kg_phase3_eur_allchrs_polymorphic.recoded.nodup.frq.gz",sep="\t",header=F,stringsAsFactors=F)
names(ref)<-c("SNP","EA","NEA","EAF")

n1<-nchar(ref$EA)
n2<-nchar(ref$NEA)

w1<-which(n1>1 & n2==1)
ref$NEA[w1]<-"D"
ref$EA[w1]<-"I"

w1<-which(n1==1 & n2>1)
ref$NEA[w1]<-"I"
ref$EA[w1]<-"D"

w1<-which(n1>1 & n2>1 & n1>n2)
ref$NEA[w1]<-"D"
ref$EA[w1]<-"I"

w1<-which(n1>1 & n2>1 & n2>n1)
ref$NEA[w1]<-"I"
ref$EA[w1]<-"D"

coh<-list()
reflist<-list()

for (i in 1:(length(l))){
#j<-as.data.frame(fread(paste(path,l[i],sep="")))
cat(i,"\n")
j<-read.table(paste(path,l[i],sep="/"),sep="\t",header=T,stringsAsFactors=F)
j<-j[order(j$CHR,j$POS,j$CpG),]
j<-j[j$CHR==20,]
assoc<-data.frame(cohort=rep(l2[i],nrow(j)),j)
assoc<-unique(assoc)
cat(nrow(assoc),"\n")

pc<-ref
m<-match(assoc$SNP,pc$SNP)
pc<-pc[na.omit(m),]
m<-match(pc$SNP,assoc$SNP)
assoc<-assoc[na.omit(m),]

cat(nrow(assoc),"\n")

#match on SNP id and effect allele to check betas
m_ea<-match(paste(pc$SNP,pc$EA,pc$NEA),paste(assoc$SNP,assoc$EA,assoc$NEA))
cat(length(m_ea),"\n")
cat(length(na.omit(m_ea)),"\n")

#match on EA vs EA
if (length(na.omit(m_ea))>0){
assoc2<-assoc[na.omit(m_ea),]
m_ea<-match(paste(assoc$SNP,assoc$EA,assoc$NEA),paste(pc$SNP,pc$EA,pc$NEA))
pc2<-pc[na.omit(m_ea),]}

#match NEA vs EA;  #swap effect allele,beta sign and EAF
m_nea<-match(paste(pc$SNP,pc$NEA,pc$EA),paste(assoc$SNP,assoc$EA,assoc$NEA))
if (length(m_nea)>0){
assoc3<-assoc[na.omit(m_nea),]
m_nea<-match(paste(assoc$SNP,assoc$EA,assoc$NEA),paste(pc$SNP,pc$NEA,pc$EA))
pc3<-pc[na.omit(m_nea),]

table(pc3$EA)

#  D   I 
#142 218

NEA<-assoc3$EA
EA<-assoc3$NEA
assoc3$EA<-EA
assoc3$NEA<-NEA
#swap beta sign
w1<-which(assoc3$BETA<0)
w2<-which(assoc3$BETA>0)

if (length(w1)>0){assoc3$BETA[w1]<-abs(assoc3$BETA[w1])}
if (length(w2)>0){assoc3$BETA[w2]<--(assoc3$BETA[w2])}
#swap EAF
assoc3$EAF[w1]<-(1-assoc3$EAF[w1])
assoc3$EAF[w2]<-(1-assoc3$EAF[w2])
}

assoc_overlap<-rbind(assoc2,assoc3) #951738
pc_overlap<-rbind(pc2,pc3)

#
#m<-match(pc_overlap$SNP,assoc_overlap$SNP)
#assoc_overlap<-assoc_overlap[m,]

str<-match(assoc$SNP,assoc_overlap$SNP)
str<-which(is.na(str))
length(str)

assoc_str<-assoc[str,]
pc_str<-pc[str,]

A<-which(assoc_str$EA=="A")
T<-which(assoc_str$EA=="T")
G<-which(assoc_str$EA=="G")
C<-which(assoc_str$EA=="C")

if (length(A)>0){assoc_str[A,"EA"]<-rep("T",length(A))}
if (length(T)>0){assoc_str[T,"EA"]<-rep("A",length(T))}
if (length(G)>0){assoc_str[G,"EA"]<-rep("C",length(G))}
if (length(C)>0) {assoc_str[C,"EA"]<-rep("G",length(C))}

A<-which(assoc_str$NEA=="A")
T<-which(assoc_str$NEA=="T")
G<-which(assoc_str$NEA=="G")
C<-which(assoc_str$NEA=="C")

if (length(A)>0){assoc_str[A,"NEA"]<-rep("T",length(A))}
if (length(T)>0){assoc_str[T,"NEA"]<-rep("A",length(T))}
if (length(G)>0){assoc_str[G,"NEA"]<-rep("C",length(G))}
if (length(C)>0) {assoc_str[C,"NEA"]<-rep("G",length(C))}

#match on SNP id and effect allele to check betas
m_ea<-match(paste(assoc_str$SNP,assoc_str$EA,assoc_str$NEA),paste(pc$SNP,pc$EA,pc$NEA))
if (length(m_ea)>0){
pc_str2<-pc[na.omit(m_ea),]

m_ea<-match(paste(pc$SNP,pc$EA,pc$NEA),paste(assoc_str$SNP,assoc_str$EA,assoc_str$NEA))
assoc_str2<-assoc_str[na.omit(m_ea),]}

m_nea<-match(paste(pc$SNP,pc$NEA,pc$EA),paste(assoc_str$SNP,assoc_str$EA,assoc_str$NEA))

if (length(m_nea)>0){
assoc_str3<-assoc_str[na.omit(m_nea),]
m_nea<-match(paste(assoc_str$SNP,assoc_str$EA,assoc_str$NEA),paste(pc$SNP,pc$NEA,pc$EA))
pc_str3<-pc[na.omit(m_nea),]

min(assoc_str3$EAF)
#[1] 0.4018
max(assoc_str3$EAF)
#[1] 0.5

#swap effect allele and beta sign
#swap effect allele
NEA<-assoc_str3$EA
EA<-assoc_str3$NEA
assoc_str3$EA<-EA
assoc_str3$NEA<-NEA

#swap beta sign
w1<-which(assoc_str3$BETA<0)
w2<-which(assoc_str3$BETA>0)
if (length(w1)>0){assoc_str3$BETA[w1]<-abs(assoc_str3$BETA[w1])}
if (length(w2)>0){assoc_str3$BETA[w2]<--(assoc_str3$BETA[w2])}

assoc_str3$EAF[w1]<-(1-assoc_str3$EAF[w1])
assoc_str3$EAF[w2]<-(1-assoc_str3$EAF[w2])

}

assoc_overlap_str<-rbind(assoc_str2,assoc_str3) #15303
pc_overlap_str<-rbind(pc_str2,pc_str3)

assoc_overlap<-rbind(assoc_overlap,assoc_overlap_str)
assoc_overlap<-unique(assoc_overlap)
pc_overlap<-rbind(pc_overlap,pc_overlap_str) #152661
pc_overlap<-unique(pc_overlap)

m<-match(assoc_overlap$SNP,pc_overlap$SNP)
pc_overlap<-pc_overlap[m,]

coh[[i]]<-assoc_overlap
reflist[[i]]<-pc_overlap

}

save(coh,reflist,file="betacomp.Robj")
#/panfs/panasas01/shared-godmc/meta-analysis/betacomparison/betacomp.Robj

