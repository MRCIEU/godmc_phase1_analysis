library(RColorBrewer)
library(ggplot2)
library(rgeos)
library(maptools)
data(wrld_simpl)

cohorts<-read.table("cohorts.txt",sep="\t",header=F)
names(cohorts)[1:2]<-c("Cohort","Country")
cohorts$value<-seq(10,400,by=20)

wrld_simpl@data$id <- wrld_simpl@data$NAME
wrld <- fortify(wrld_simpl, region="id")
wrld <- subset(wrld, id != "Antarctica") # we don't rly need Antarctica


getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cohorts$cols<-getPalette(length(cohorts$Country))
o<-order(getPalette(length(cohorts$Country)),decreasing=T)
cohorts$cols<-cohorts$cols[o]
cohorts$Country<-as.factor(cohorts$Country)
o<-order(cohorts$Country)
cohorts<-cohorts[o,]

w<-which(cohorts$Country%in%c("United Kingdom"))

spl<-strsplit(as.character(cohorts$Cohort[w]),split=",")
l<-length(spl[[1]])/2
uk1<-paste(spl[[1]][1:l],sep="",collapse=", ")
uk2<-paste(spl[[1]][(l+1):(length(spl[[1]]))],collapse=", ",sep="")

cohorts$Cohort<-as.character(cohorts$Cohort)
addrow<-cohorts[w,]
addrow$Cohort<-uk2
addrow$Country<-"United Kingdom2"
cohorts<-rbind(cohorts,addrow)
cohorts[w,"Cohort"]<-uk1
o<-order(as.character(cohorts$Country))
cohorts<-cohorts[o,]

gg <- ggplot()
gg <- gg + geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="white", color="#7f7f7f", size=0.4)
gg <- gg + geom_map(data=cohorts, map=wrld, aes(map_id=Country,fill=Country),  color="white", size=0.4)
 # this sets the scale and, hence, the legend
gg <- gg + scale_fill_manual(values=cohorts$cols,name="Cohort",labels=cohorts$Cohort)
#gg<-gg+scale_fill_brewer(palette="Set1",name="Cohort",labels=cohorts$Cohort)
gg<-gg + theme(legend.position="bottom")
gg<-gg+guides(fill=guide_legend(nrow=11,title.position="top",size=0.5))
gg<-gg+theme(legend.text = element_text(size = 6))
gg<-gg+ theme(legend.title = element_text(size = 8))
gg<-gg+scale_x_discrete(breaks=NULL,name="")
gg<-gg+scale_y_discrete(breaks=NULL,name="")
ggsave(gg, file="map.pdf",width=8,height=6)
