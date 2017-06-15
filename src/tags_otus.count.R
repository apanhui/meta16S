#!/usr/bin/Rscript

# =======================
# linpeng
# plin@genedenovo.com
# Wed Jun 14 11:16:34 CST 2017
# =======================

args<-commandArgs(T)
if(length(args)!=2){
	cat("Rscript two_axes.R <file> <opre>\n")
	q()
}
data<-read.table(args[1],sep="\t",header=T)
.libPaths("/home/linpeng/R/x86_64-unknown-linux-gnu-library/3.2")
library(reshape2)
library(RColorBrewer)
library(TeachingDemos)

##data1 was the barplot data
data1<-data[data[,1]!="Avg",]
#convert data1 to matrix
data1<-dcast(data1,SampleID~attr)
row.names(data1)<-data1[,1]
bar_order<-c("Total_Tags","Taxon_Tags","Unclassified_Tags",
             "Unique_Tags","Singleton_Tags","OTUs")
data1<-data1[,bar_order]
head(data1)

##data2 was the legend data
data2<-data[data[,1]=="Avg",]
row.names(data2)<-data2$attr
data2<-data2[bar_order,]
legend_text<-apply(data2,1,
                   function(x){paste(x[2],"(",as.integer(x[3]),")",sep="")})

##set the axes's range and change Otu num
y1<-ceiling(max(data1[,"Total_Tags"])/10000)*10000+10000
y2<-ceiling(max(data1[,"OTUs"]/100))*100+100
data3<-data1
data3$OTUs<-data1$OTUs*(y1/y2)

##calculate pictures' width
if(nrow(data1)<=14){
      pic_wd<-12
}else{
      pic_wd<-(nrow(data1)-14)*0.3+12
}
cat("pictures' width :",pic_wd,"inch\n")

opar<-par(no.readonly = T)

##color 
#Spectral<-brewer.pal(6,"Spectral")
Set2<-mypalette<-brewer.pal(6,"Set2")
#topo<-c("#F2414A","#B32E90","#FFF540","#2BC482","#F5AF3F","2D91E3")
Spectral=c("#F2414A","#B32E90","#FFF540","#2BC482","#F5AF3F","#2D91E3")

### png plot ###
png(file=paste(args[2],"png",sep="."), width=pic_wd,height=8,units="in",res=300)
par(mar=c(5,4,4,4)+0.1)
bar1<-barplot(t(data3),beside=T,col=Spectral,ylim=c(0,y1),
        axes=T,xlab="Sample Names",ylab="Number of Tags")
legend("topleft",inset=c(0,0),legend_text,fill=Spectral,bty='n',
      x.intersp=0.1,cex=1,xpd=T,ncol=3) #delete text.width=5
##plot number over the bar
par(srt=90)
text(x=bar1,y=as.vector(t(data3))+200,lab=as.vector(t(data1)),adj=0,cex=0.6)
##set the second axis and lable
updateusr(1:2,range(0,y1),1:2,range(0,y2))
axis(4)
mtext("Number of OTUs",side=4,padj=5)
dev.off()

### pdf plot ###
pdf(file=paste(args[2],"pdf",sep="."), width=pic_wd,height=8)
par(mar=c(5,4,4,4)+0.1)
bar1<-barplot(t(data3),beside=T,col=Spectral,ylim=c(0,y1),
        axes=T,xlab="Sample Names",ylab="Number of Tags")
legend("topleft",inset=c(0,0),legend_text,fill=Spectral,bty='n',
      x.intersp=0.1,cex=1,xpd=T,ncol=3) #delete text.width=8
##plot number over the bar
par(srt=90)
text(x=bar1,y=as.vector(t(data3))+200,lab=as.vector(t(data1)),adj=0,cex=0.6)
##set the second axis and lable
updateusr(1:2,range(0,y1),1:2,range(0,y2))
axis(4)
mtext("Number of OTUs",side=4,padj=4)
dev.off()
