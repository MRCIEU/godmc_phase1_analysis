library(meffil)
args <- (commandArgs(TRUE));
user <- as.character(args[1]);
cohort <- as.character(args[2]);


y<-meffil.featureset("450k")

path="/panfs/panasas01/shared-godmc/clump"

coh<-read.table("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/extract_sftp/data/cohorts.txt",sep=" ",header=F)

coh.out<-data.frame(CpG=y$name,CHR=y$chromosome,POS=y$position)
y<-coh.out

getcistrans <- function(coh.out, res, win)
{
	a <- subset(coh.out, CpG==res[1])
	b <- do.call(rbind, strsplit(as.character(res[-1]), split=":"))
	b <- data.frame(b, stringsAsFactors=FALSE)
	b[,2] <- as.numeric(b[,2])
	b$cis <- b[,1] == a$CHR & b[,2] > a$POS - win & b[,2] < a$POS + win
	return(b)
}

#getcistrans(coh.out, spl[[1]], 1000000)

p<-paste(path,"/",user,"_",cohort,"/",sep="")
l<-list.files(path=p)
g<-grep("indexSNP",l)
l<-l[g]

r.out<-data.frame()

for (k in 1:length(l)){ 
r<-readLines(paste(p,l[k],sep=""))
cat(l[k],"\n")
spl<-strsplit(r,split=" ")


cistrans<-lapply(spl,function(obj) x<-getcistrans(coh.out, obj, 1000000))
cis<-lapply(cistrans,function(obj) x<-length(which(obj$cis==TRUE)))
trans<-lapply(cistrans,function(obj) x<-length(which(obj$cis==FALSE)))
cis<-do.call("rbind",cis)
trans<-do.call("rbind",trans)

snpno<-lapply(spl,function(obj) x<-(length(obj)-1))
snpno<-do.call("rbind",snpno)
probe<-lapply(spl,function(obj) x<-obj[1])
probe<-do.call("rbind",probe)

r2<-data.frame(CpG=probe,snpno=snpno, cis=cis,trans=trans)
r.out<-rbind(r.out,r2)
}

f<-paste(user,"_",cohort,".numberofindependentloci.Robj",sep="")

save(r.out,file=f)



