#!/Bio/bin/Rscript
#-------------------------------------------------------------
# name:    pca.R
# func:    do PCA analysis and draw PCA 2D & 3D plots
# version: 1.0
# author:  aipeng (pai@genedenovo.com)
#-------------------------------------------------------------
# import packages
library(gmodels)
library(ggplot2)
library(scatterplot3d)

args<-commandArgs(T)

if (length(args)  < 1)
{
  usage = "
Version: 1.0 

Descs: An R script to do PCA analysis and draw 2D & 3D scatter plots 

Usage: Rscript pca.R <matrix file> [start,end] [group.info]

Note: 1. the matrix file must be seperated by <TAB>
      2. the first line of the matrix file is the sample names 
      3. the first column of the matrix file is the gene/tag names 
      4. [strat,end] is used to define the data start and end column, default will use the whole matrix 
      5. [group.info] was defined like : A,B,B,B,C,C, seperated by comma, relative to sample names one by one
\n"
  cat(usage)
  quit("no")
}

# read the data
data = read.table(args[1],header=T,row.names=1,sep="\t",check.names=F)

# cut the target data if defined
if(length(args)>=2)
{
  colms<-strsplit(args[2],",")

  start=as.numeric(colms[[1]][1])-1
  end=as.numeric(colms[[1]][2])-1

  data=data[,start:end]
  data=data[which(rowSums(data)>0),]

  cat("[LOG] fetch the specific columns of the matrix:",args[2],"\n")
  cat("[LOG] sum rows:",nrow(data),"\n")
  cat("[LOG] sum samples:",ncol(data),"\n")
}

# the samples less than 2 will exit the program 
if(nrow(data)<2)
{
  cat("\n###\nthe taxa or otu count less than 2, can not draw pca 2D,quit\n###\n")
  quit("no")
}

# prepare the data for PCA analysis
# do scale for the data 
tmp <- matrix(0,nrow = nrow(data), ncol = ncol(data))
for(i in 1:ncol(data)){
  tmp[,i] = data[ ,i]/sum(data[ ,i])
}
# tmp<-t(scale(t(tmp)))
colnames(tmp) <- colnames(data)
rownames(tmp) <- rownames(data)
data <- t((tmp))

#----------------------------------------------------
# do PCA analysis 
data.pca <- fast.prcomp(data,retx=T,scale=F,center=T)

cat("[LOG] PCA analysis was done ...\n")

# parse the PCA analysis result 
summary <- summary(data.pca)
result <- summary$importance

# save the PCA result 
rotation = data.pca$rotation
write.table(cbind(ID=rownames(rotation),rotation),file=paste(args[1],".PCA.component.xls",sep=""),sep="\t", quote=FALSE,row.names=FALSE)

x = data.pca$x
write.table(cbind(Samples=rownames(x),x),file=paste(args[1],".PCA.samples.xls",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

print(max(x[,1]))

cat("[LOG] the result of PCA was saved into files:\n\t",
  paste(args[1],"PCA.samples.xls",sep=""),"and\n\t",paste(args[1],".PCA.component.xls",sep=""),"\n")

# fetch the PC value 
pro1 = as.numeric(sprintf("%.3f",result[2,1]))*100
pro2 = as.numeric(sprintf("%.3f",result[2,2]))*100
pro3 = 0

if(length(result[2,])>2)
{
  pro3 = as.numeric(sprintf("%.3f",result[2,3]))*100  
}

# sample names 
sample_names = rownames(data)
cat("[LOG] sample names:",sample_names,"\n")

# create the data for 2D PCA plot 
# turn matrix to frame 
pc = as.data.frame(summary$x)
pc$names = sample_names

# read the group information
group = sample_names
legend_title = "Samples"

# creat the color for 3D plot 
# the function to create default ggplot2 colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

sample_colour = gg_color_hue(length(sample_names))

if (length(args) > 2)
{
  group = unlist(strsplit(args[3],","))
  legend_title = "Group"
  cat("[LOG] group information was read:",group,"\n")

  group_color = gg_color_hue(length(unique(group)))
  group_frame = data.frame(color=group_color)
  rownames(group_frame) = unique(group)
  
  for (i in 1:length(sample_names))
  {
    group_name = group[i]
    sample_colour[i] = as.character(group_frame[group_name,1])
  }
}

pc$group = group

#-------------------------------------------------------------
# draw 2D PCA with ggplot2
xlab=paste("PC1(",pro1,"%)",sep="") 
ylab=paste("PC2(",pro2,"%)",sep="")

pca = ""
noname_pca = ""

if (length(unique(pc$group)) > 6)
{
pca=ggplot(pc,aes(PC1,PC2)) + 
  geom_point(size=3,aes(color=group)) + 
  geom_text(aes(label=names),size=4,vjust=-1) +
  labs(x=xlab,y=ylab,title="PCA",color=legend_title,shape=legend_title) + 
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + 
  theme_bw()

nonames_pca=ggplot(pc,aes(PC1,PC2)) + 
  geom_point(size=3,aes(color=group)) + 
  labs(x=xlab,y=ylab,title="PCA",color=legend_title,shape=legend_title) + 
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + 
  theme_bw()
} else {
pca=ggplot(pc,aes(PC1,PC2)) + 
  geom_point(size=3,aes(shape=group,color=group)) + 
  geom_text(aes(label=names),size=4,vjust=-1) +
  labs(x=xlab,y=ylab,title="PCA",color=legend_title,shape=legend_title) + 
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + 
  theme_bw()

nonames_pca=ggplot(pc,aes(PC1,PC2)) + 
  geom_point(size=3,aes(shape=group,color=group)) + 
  labs(x=xlab,y=ylab,title="PCA",color=legend_title,shape=legend_title) + 
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + 
  theme_bw()
}

ggsave(paste(args[1],".PCA2D.pdf",sep=""),pca,width=10,height=8)
ggsave(paste(args[1],".PCA2D.png",sep=""),pca,width=10,height=8)

cat("[LOG] 2D PCA plot with sample names was created :)\n")

ggsave(paste(args[1],".PCA2D.nonames.pdf",sep=""),nonames_pca,width=10,height=8)
ggsave(paste(args[1],".PCA2D.nonames.png",sep=""),nonames_pca,width=10,height=8)

cat("[LOG] 2D PCA plot without sample names was created :)\n")

#-------------------------------------------------------------
# draw 3D PCA with scatterplot3d
# fetch the PC limits 
xmax = max(data.pca$x[,1])
xmin = min(data.pca$x[,1])
ymax = max(data.pca$x[,2])
ymin = min(data.pca$x[,2])
zmin = min(data.pca$x[,3])
zmax = max(data.pca$x[,3])

if(abs(zmax-zmin)<=0.0005)
{
  zmax=0.001
  zmin=-0.001
}

zlab = paste("PC3(",pro3,"%)",sep="")

legend_ncol = as.integer((length(unique(sample_colour)) - 1)/30) + 1
inset = -0.2*legend_ncol + 0.15
right_margin = 2 + 8*legend_ncol

# A function for draw 3D PCA plots 
draw_3d_pca = function(file,width=10,height=8,names=TRUE,legend=TRUE) {
  suffix = unlist(strsplit(file,".",fixed=TRUE))
  suffix = suffix[length(suffix)]
  
  # open dev 
  if (suffix == "pdf") {
    pdf(file,width=width,height=height)
  } else {
    png(file,width=width,height=height)
  }
  
  # draw 3D plot 
  s3d = scatterplot3d(data.pca$x[,1],data.pca$x[,2],data.pca$x[,3],
    main="PCA",xlab=xlab,ylab="",zlab=zlab,
    mar=c(5,3,4,right_margin)+0.1,
    xlim=c(1.2*xmin,1.2*xmax),ylim=c(1.2*ymin,1.2*ymax),zlim=c(1.2*zmin,1.2*zmax),
    color=sample_colour,angle=35,pch=16,type="h",scale.y=0.8,cex.symbols=1.2)
  
  # add the y lab by text (the default y lab position is terrible)
  dims <- par("usr")
  x <- dims[1]+ 0.9*diff(dims[1:2])
  y <- dims[3]+ 0.08*diff(dims[3:4])
  text(x,y,ylab,srt=45)

  ## add text to 3D plot 
  if (names == TRUE) {
    text(s3d$xyz.convert(data.pca$x[,1],data.pca$x[,2],data.pca$x[,3]),labels=sample_names,pos=3,cex=1.2)
  }
  
  ## add the legend for 3D plot 
  if (length(args) > 2 && legend == TRUE)
  {
    group = unlist(strsplit(args[3],","))
    group_names = unique(group)
    group_color = gg_color_hue(length(group_names))
    
    legend("right",inset=inset,group_names,col=group_color,pch=16,cex=1.2,title=legend_title,bty="n",xpd=TRUE,ncol=legend_ncol)
  } else if (legend == TRUE) {
    legend("right",inset=inset,sample_names,col=sample_colour,pch=16,cex=1.2,title=legend_title,bty="n",xpd=TRUE,ncol=legend_ncol)
  }

  # close dev
  dev.off()
}

width = 10 + legend_ncol

#-------------------------------------
# draw 3D plot with sample names (pdf)
draw_3d_pca(file=paste(args[1],".PCA3D.pdf",sep=""),names=TRUE,legend=TRUE,width=width,height=8)

#-------------------------------------
# draw 3D plot without sample names (pdf)
draw_3d_pca(file=paste(args[1],".PCA3D.nonames.pdf",sep=""),names=FALSE,legend=TRUE,width=width,height=8)

cat("[LOG] 3D PCA plot (pdf) was created :)\n")

#-------------------------------------
# draw 3D plot with sample names (png)
draw_3d_pca(file=paste(args[1],".PCA3D.png",sep=""),names=TRUE,legend=TRUE,width=width*100,height=800)

#-------------------------------------
# draw 3D plot without sample names (png)
draw_3d_pca(file=paste(args[1],".PCA3D.nonames.png",sep=""),names=FALSE,legend=TRUE,width=width*100,height=800)

cat("[LOG] 3D PCA plot (png) was created :)\n")

# end of the program
quit("no")
