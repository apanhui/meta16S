#!/usr/bin/Rscript

# =======================
# linpeng
# plin@genedenovo.com
# Mon Apr 24 11:07:48 CST 2017
# =======================

args<-commandArgs(T)
if(length(args)<2){
	cat("Rscript heatmap.R <file> <out_prefix>\n")
	q()
}
library(pheatmap)
library(gplots)

data<-read.table(args[1],head = TRUE,sep="\t",row.names=1,check.names=F)
data1=data

if(nrow(data1)<2)
{
	cat("\n###\nthe taxa or otu count less than 2, can not draw heatmap,quit\n###\n")
	quit("no")
}

data1=t(scale(t(data1)))
fontsize_row=10
fontsize_col=10

max_row_len=0.6*fontsize_row*max(nchar(rownames(data1)))
width_pt=70+(10*ncol(data1))+max_row_len
if(width_pt>800){
  width_pt<-width_pt
}else if(width_pt<600){
  cell_wpt<-(width_pt-70-max_row_len)/ncol(data1)
#  cat("cell_wpt :",cell_wpt,"\n")
  while(fontsize_col<25){
      fontsize_col<-fontsize_col+1
#      cat("fontsize_col :",fontsize_col,"\n")
    if( (cell_wpt/2)<(fontsize_col*0.8) ){
      break
    }
  }
  if(cell_wpt>40){
    width_pt<-70+(40*ncol(data1))+max_row_len
  }
}else{
  width_pt<-800
}
width_pt<-ceiling(width_pt/10)*10
cat("width_pt :",width_pt,"\n")

max_col_len=0.6*fontsize_col*max(nchar(colnames(data1)))
height_pt=50+(10*nrow(data1))+max_col_len
if(height_pt>1200){
  while(fontsize_row>4){
    fontsize_row<-fontsize_row-1;
    height_pt<-50+(fontsize_row*nrow(data1))+max_col_len
    if(height_pt<1200){
      height_pt<-1200
      break
    }
  }
}else if(height_pt<800){
  height_pt=height_pt
}else{
  height_pt<-1200
}
height_pt<-ceiling(height_pt/10)*10
cat("height_pt :",height_pt,"\n")
cat("fontsize_row :",fontsize_row,"\n")
cat("fontsize_col :",fontsize_col,"\n")

pdf(paste(args[2],".heatmap.pdf",sep=""),ceiling(width_pt/75),ceiling(height_pt/75))
#heatmap.2(as.matrix(data1),col=greenred(255),revC=FALSE, scale="row",key=TRUE,symkey=TRUE,trace="none", density.info="none",cexRow=0.7,cexCol=0.8,keysize=0.8,margins=c(8,8))
pheatmap(as.matrix(data1), scale = "row",fontsize_row=fontsize_row,fontsize_col=fontsize_col)
dev.off()
cat("output file : ",args[2],".heatmap.pdf","\n",sep="")

png(paste(args[2],".heatmap.png",sep=""),width_pt,height_pt)
#heatmap.2(as.matrix(data1),col=greenred(255),revC=FALSE, scale="row",key=TRUE,symkey=TRUE,trace="none", density.info="none",cexRow=0.7,cexCol=0.8,keysize=0.8,margins=c(8,8))
pheatmap(as.matrix(data1), scale = "row",fontsize_row=fontsize_row,fontsize_col=fontsize_col)
dev.off()
cat("output file : ",args[2],".heatmap.png","\n",sep="")

if(fontsize_row==4){
	height_pt=ceiling((50+(2*nrow(data1))+max_col_len)/10)*10
	width_pt=ceiling((width_pt-max_row_len)/10)*10
	pdf(paste(args[2],".without_rowname.heatmap.pdf",sep=""),ceiling(width_pt/75),ceiling(height_pt/75))
	pheatmap(as.matrix(data1), scale = "row",fontsize_col=fontsize_col,show_rownames=F)
	dev.off()
	cat("output file : ",args[2],".without_rowname.heatmap.pdf","\n",sep="")
	png(paste(args[2],".without_rowname.heatmap.png",sep=""),width_pt,height_pt)
	pheatmap(as.matrix(data1), scale = "row",fontsize_col=fontsize_col,show_rownames=F)
	dev.off()
	cat("output file : ",args[2],".without_rowname.heatmap.png","\n",sep="")
}
quit("no")
