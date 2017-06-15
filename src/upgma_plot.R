#!/Bio/bin/Rscript
#-----------------------------------------------------
# name: upgma_plot.R
# func: draw UPGMA tree and taxonomy stack fill bar
# version: 1.0
# date: 2017-05-27
#-----------------------------------------------------
args<-commandArgs(T)

if (length(args) < 5){
	usage = "
Version: 1.0 

Usage:   Rscript upgma_plot.R <tree file> <taxonomy abudance file> <Taxonomy Level> <Weighted|UnWeighted> <size> [groups info]
\n"
	cat(usage)
	quit("no")
}

# import packages 
library(vegan)
library(ape)

# init stack bar color 
col=c("red","green3","blue","cyan","magenta","yellow","gray","purple","brown",
      "darkred","darkmagenta","mediumvioletred","peru","lightcoral","darksalmon",
      "sandybrown","black","darkgreen","forestgreen","mediumseagreen","limegreen",
      "darkkhaki","mediumspringgreen","lawngreen","lightgreen","greenyellow","navy",
      "midnightblue","darkviolet","darkslateblue","blueviolet","dimgray","lightslateblue",
      "lightslategray","cornflowerblue","lightseagreen","violet","darkgray",
      "darkturquoise","mediumturquoise","lightgray","powderblue","gainsboro",
      "palegoldenrod","moccasin","lavender","blanchedalmond","papayawhip","linen",
      "beige","oldlace","aliceblue","lightgoldenrodyellow","ghostwhite","floralwhite",
      "mintcream")

tree_file = args[1]
taxon_abu = args[2]
outname   = args[3]
distance  = args[4]
size      = as.numeric(args[5])

distance_name = paste(distance,"Unifrac Distance",sep=" ")

# prepare UPGMA tree data
mytree <- read.tree(tree_file)

# prepare stack fill bar taxonomy data
taxon_data = read.table(taxon_abu,header=TRUE,sep="\t",check.names=F,comment.char = "")
rownames(taxon_data) = taxon_data[,1]
taxon_data = taxon_data[,-1]

# create the stack data by the sample order of the tree 
samples = mytree$tip.label
taxones = rownames(taxon_data)
samples_count = length(samples)

stack_data = matrix(0,nrow = nrow(taxon_data) , ncol = ncol(taxon_data))
rownames(stack_data) = taxones
colnames(stack_data) = samples

taxon_data[,samples[2]]
for(i in 1:ncol(taxon_data)){
  sample = samples[i]
  stack_data[,i] = taxon_data[ , sample] / sum( taxon_data[ ,sample] )
}

# set the edge and sample colors
edge_color   = rep(1,nrow(mytree$edge))
sample_color = rep(1,samples_count)

mark=1     ## dont draw branch color

if(length(args)>=6)
{
  colms<-strsplit(args[6],",")
  classify<-as.numeric(colms[[1]])
  k=1
  j=1
  if(mark!=1)
  {
    for(i in 1:nrow(mytree$edge))
    {
      if(mytree$edge[i,2]!=classify[j])
      {
        edge_color[i]=col[k]
      }else{
        edge_color[i]=col[k]
        k=k+1
        j=j+1
      }
    }
  }
  
  cnt<-1
  for(i in 1:length(classify))
  {
    for(j in 1:classify[i] )
    {
      sample_color[cnt]<-col[i]
      cnt<-cnt+1
    }
  }
}

mat = t(matrix(c(1,2),nrow=2))
xrange=node.depth.edgelength(mytree)[1]

# draw pdf figure
pdf_width  = 16
pdf_height = 8

pdf(paste(distance,outname,"UPGMA_stack.pdf",sep="."),pdf_width,pdf_height)
layout(mat,widths=c(2,3))
par(mar=c(7,3,3,0.8))

plot.phylo(mytree,label.offset=1,edge.width=2,edge.color=edge_color,tip.color=sample_color,
           y.lim=c(1-0.45,samples_count+0.45),show.tip.label=F)
text(x=xrange*1.05,y=1:samples_count,labels=samples,xpd=TRUE)

# add the tree distance axis 
axis(1,at=seq(0,xrange,by=0.05),labels=seq(0,xrange,by=0.05),lwd=0,lwd.ticks=1)
lines(c(0,xrange),c(0.28,0.28),xpd=T)

text(x=xrange/2,y=0.55-(samples_count)/7,labels=distance_name,xpd=T)

# draw stack part
par(mar=c(7,0.8,3,size))
dd=barplot(as.matrix(stack_data),xlab="Relative Abundance",main="",xlim=c(0,1),
           col=col[1:nrow(stack_data)],las=1,horiz=T,axisnames=F)

legend(x=1.03,y=median(dd),yjust=0.5,legend=taxones,cex=1.1,pt.cex=2.5,
       col=col[1:nrow(stack_data)],pch=15,xpd=T,title=outname)

dev.off()

# draw png figure
png(paste(distance,outname,"UPGMA_stack.png",sep="."),width=1600,height=800)
layout(mat,widths=c(2,3,lcm(size)))

par(mar=c(7,3,3,0.8))
plot.phylo(mytree,label.offset=1,edge.width=2,edge.color=edge_color,tip.color=sample_color,
           y.lim=c(1-0.45,samples_count+0.45),show.tip.label=F)
text(x=xrange*1.05,y=1:samples_count,labels=samples,xpd=TRUE)

# add the tree distance axis 
axis(1,at=seq(0,xrange,by=0.05),labels=seq(0,xrange,by=0.05),lwd=0,lwd.ticks=1)
lines(c(0,xrange),c(0.28,0.28),xpd=T)

text(x=xrange/2,y=0.55-(samples_count)/7,labels=distance_name,xpd=T)

# draw stack part
par(mar=c(7,0.8,3,size))
dd=barplot(as.matrix(stack_data),xlab="Relative Abundance",main="",xlim=c(0,1),
           col=col[1:nrow(stack_data)],las=1,horiz=T,axisnames=F)

legend(x=1.03,y=median(dd),yjust=0.5,legend=taxones,cex=1.1,pt.cex=2.5,
       col=col[1:nrow(stack_data)],pch=15,xpd=T,title=outname)

dev.off()

