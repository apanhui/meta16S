#!/Bio/bin/Rscript
#-------------------------------------------------------------
# name:    otus_venn.R
# func:    draw otus venn diagram of different groups (batch mode)
# version: 1.0
# author:  aipeng (pai@genedenovo.com)
#-------------------------------------------------------------
args<-commandArgs(T)
if (length(args) != 2){
  usage = "
Version: 1.0

Descs:   An R script to batch draw venn diagram of compare different group OTUs

Usage:   Rscript otus_venn.R <all.otus.group.shared> <compare str>

Note:    compare str was defined the in the meta16 conf file, (two_group_diff or multi_group_diff)
\n"
  cat(usage)
  quit("no")
}

# import packages
library(VennDiagram)
library(UpSetR)

# int the fill colors 
my_colors = c("cornflowerblue","green","darkorchid1","purple","yellow")

# read the group shared info
data = read.table(file=args[1],header = TRUE,sep="\t")
data = t(data[,-1])
data = data[-2,]
colnames(data) = data[1,]
data = data[-1,]

# all otus and  all groups 
otus = rownames(data)
groups = colnames(data)

# read the compare
diffs = args[2]
pairwises = unlist(strsplit(diffs,","))

# a function to draw venn diagram, create png (600 dpi) and pdf two format file
venn_plot = function (x,outname,cex=1.5){
  fill_colors = my_colors[1:length(x)]
  
  # png format 
  filename = paste(outname,".venn.png",sep="")
  venn.diagram(x,filename=filename,col="white",fill=fill_colors,lwd=.5,
                           cex=cex,alpha=.5,resolution=600,width=1200,height=1200,imagetype="png")
  # pdf format
  venn = venn.diagram(x,filename=NULL,col="white",fill=fill_colors,
                           cex=cex,alpha=.5)
  pdf_name = paste(outname,".venn.pdf",sep="")
  pdf(file=pdf_name,height=7,width=7)
  grid.draw(venn)
  dev.off()
}

# a function to draw venn diagram with UpSetR and output the elements of each sets
upset_plot = function (x,outname){
  x.all = sort(unique(unlist(x)))
  x.mat = matrix(0, nrow=length(x.all), ncol=length(x), 
    dimnames=list(genes=x.all, groups=names(x)))
  for(i in 1:length(x)) x.mat[x[[i]], i] <- 1
  
  png(paste(outname,".upset.png",sep=""),width=2000,height=1500,res=300)
  upset(as.data.frame(x.mat),nsets=length(x),nintersects=NA)
  dev.off

  pdf(paste(outname,".upset.pdf",sep=""),width=8,height=6)
  upset(as.data.frame(x.mat),nsets=length(x),nintersects=NA)
  dev.off

  # output the elements
  x.venn <- split(rownames(x.mat), apply(x.mat, 1, paste, collapse=""))
  names(x.venn) <- 
    sapply(names(x.venn), function(.ele) 
    paste(colnames(x.mat)[as.logical(as.numeric(strsplit(.ele, "")[[1]]))], 
    collapse="&"))
  
  names    = names(x.venn)
  numbers  = rep(0,length(names))
  elements = rep("",length(names))

  for (i in 1:length(names)){
    numbers[i] = length(x.venn[[i]])
    elements[i] = paste(x.venn[[i]],collapse=",")
  }

  outdata = data.frame(Group=names,Number=numbers,Elements=elements)
  write.table(outdata,file=paste(outname,".venn.elements.xls",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
}

for (i in pairwises){
  compares = unlist(strsplit(i,"&"))
  x = list()
  for (gname in compares){
    x[[gname]] =otus[as.numeric(data[,gname]) > 0]
  }
  name = gsub(pattern="&",replacement="_vs_",i)
  
  if (length(x) <= 5){
    venn_plot(x,name)
  } else {
    warning("VennDiagram can't draw >5 sets venn diagram, skip ...")
  }

  upset_plot(x,name)
}
