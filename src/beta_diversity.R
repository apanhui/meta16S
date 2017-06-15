#!/Bio/bin/Rscript
#-------------------------------------------------------------
# name:    beta_diversity.R
# func:    1. calc (un)weighted unifrac distance
#          2. do PCoA
#          3. do NMDS
# version: 1.0
# author:  aipeng (pai@genedenovo.com)
# date:    2017-05-24
#-------------------------------------------------------------
args = commandArgs(T)

if (length(args) != 3)
{
  cat("
Version: 1.0

Descs: An R script to calc weighted and unweighted unifrac distance and do some relative analysis.
       1. do PCoA analysis and draw 2D scatter plot 
       2. do NMDS analysis and draw 2D scatter plot 

Usage: Rscript beta_diversity.R <otus.abundance.tsv> <otus.tree> <group file>
\n")
  q()
}

# import package 
library(GUniFrac)
library(ggplot2)
library(pheatmap)

expr_file = args[1]
tree_file = args[2]
group_file = args[3]

group_data = read.table(file=group_file, sep="\t", header=T, check.names=F, comment.char = "")
groups = group_data$Group

expr = read.table(file=expr_file, row.names=1, sep="\t", header=T, check.names = F)
tree = read.tree(file=tree_file)

#-------------------------------------------------------------
# step1: 
# calc unifrac
otu.tab.rff <- Rarefy(t(expr))$otu.tab.rff
unifracs <- GUniFrac(otu.tab.rff,tree,alpha=c(0, 0.5, 1))$unifracs

# Weighted UniFrac
wu<-unifracs[,,"d_1"]

# Unweighted UniFrac
uw<-unifracs[,,"d_UW"]

wu.out<-cbind(row.names(wu),wu);
uw.out<-cbind(row.names(uw),uw);

# save unifrac result to tsv file
write.table(wu.out,file = "weighted_unifrac.tsv",sep ="\t",quote=F,row.names=F)
write.table(uw.out,file = "unweighted_unifrac.tsv",sep ="\t",quote=F,row.names=F)

#-------------------------------------------------------------
# step2:
# draw heatmap 
pheatmap(wu,scale="none",cluster_cols=FALSE,cluster_rows=FALSE,
         cellwidth=40,cellheight=40,display_numbers=TRUE,number_format="%.4f",
         main="weighted unifrac heatmap",filename="weighted_unifrac.heatmap.png")
pheatmap(wu,scale="none",cluster_cols=FALSE,cluster_rows=FALSE,
         cellwidth=40,cellheight=40,display_numbers=TRUE,number_format="%.4f",
         main="weighted unifrac heatmap",filename="weighted_unifrac.heatmap.pdf")

pheatmap(uw,scale="none",cluster_cols=FALSE,cluster_rows=FALSE,
         cellwidth=40,cellheight=40,display_numbers=TRUE,number_format="%.4f",
         main="unweighted unifrac heatmap",filename="unweighted_unifrac.heatmap.png")
pheatmap(uw,scale="none",cluster_cols=FALSE,cluster_rows=FALSE,
         cellwidth=40,cellheight=40,display_numbers=TRUE,number_format="%.4f",
         main="unweighted unifrac heatmap",filename="unweighted_unifrac.heatmap.pdf")

#-------------------------------------------------------------
# step3:
# do PCoA
pcoa_plot = function(x,outname,k=2,title=NULL){
  pcoa = cmdscale(x,k=k,eig=T)
  
  eig = pcoa$eig
  pco1 = round(eig[1]/sum(eig)*100,2)
  pco2 = round(eig[2]/sum(eig)*100,2)
  
  pcoa_points = as.data.frame(pcoa$points)
  
  xlab = paste("PCO1 (",pco1,"%)",sep="")
  ylab = paste("PCO2 (",pco2,"%)",sep="")
  
  pcoa.plot = ""
  pcoa.plot.nonames = ""

  if (length(unique(groups)) > 6){
    pcoa.plot = ggplot(pcoa_points,aes(V1,V2)) + 
      geom_point(size=3,aes(color=groups)) + 
      geom_text(aes(label=rownames(pcoa_points)),size=4,vjust=-1)+
      labs(x=xlab,y=ylab,title=title) + 
      geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
  
    pcoa.plot.nonames = ggplot(pcoa_points,aes(V1,V2)) + 
      geom_point(size=3,aes(color=groups)) + 
      labs(x=xlab,y=ylab,title=title) + 
      geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
  } else { 
    pcoa.plot = ggplot(pcoa_points,aes(V1,V2)) + 
      geom_point(size=3,aes(shape=groups,color=groups)) + 
      geom_text(aes(label=rownames(pcoa_points)),size=4,vjust=-1)+
      labs(x=xlab,y=ylab,title=title) + 
      geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
  
    pcoa.plot.nonames = ggplot(pcoa_points,aes(V1,V2)) + 
      geom_point(size=3,aes(shape=groups,color=groups)) + 
      labs(x=xlab,y=ylab,title=title) + 
      geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
  }

  ggsave(paste(outname,"PCoA.png",sep="."),pcoa.plot,width=10,height=8)
  ggsave(paste(outname,"PCoA.pdf",sep="."),pcoa.plot,width=10,height=8)
  
  ggsave(paste(outname,"PCoA.nonames.png",sep="."),pcoa.plot.nonames,width=10,height=8)
  ggsave(paste(outname,"PCoA.nonames.pdf",sep="."),pcoa.plot.nonames,width=10,height=8)
}

pcoa_plot(wu,"weighted_unifrac",title="Weighted unifrac PCoA")
pcoa_plot(uw,"unweighted_unifrac",title="Unweighted unifrac PCoA")

#-------------------------------------------------------------
# step4:
# do NMDS

nmds_plot = function(x,outname,k=2,title=NULL){
  nmds = metaMDS(x)
  nmds.points = as.data.frame(nmds$points)
  stress<-paste("stress =",round(nmds$stress,3))

  nmds.plot = ""
  nmds.plot.nonmaes = ""

  if (length(unique(groups)) > 6){
    nmds.plot<-ggplot(nmds.points,aes(MDS1,MDS2)) + 
      geom_point(size=3,aes(color=groups)) + 
      geom_text(aes(label=rownames(nmds.points)),size=4,vjust=-1)+ 
      labs(x="NMDS1",y="NMDS2",title=title) +
      annotate("text",label=stress,x=max(nmds.points$MDS1)*.9,y=max(nmds.points$MDS2),size=5)+
      geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
  
    nmds.plot.nonmaes<-ggplot(nmds.points,aes(MDS1,MDS2)) + 
      geom_point(size=3,aes(color=groups)) + 
      labs(x="NMDS1",y="NMDS2",title=title) +
      annotate("text",label=stress,x=max(nmds.points$MDS1)*.9,y=max(nmds.points$MDS2),size=5)+
      geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    nmds.plot<-ggplot(nmds.points,aes(MDS1,MDS2)) + 
      geom_point(size=3,aes(shape=groups,color=groups)) + 
      geom_text(aes(label=rownames(nmds.points)),size=4,vjust=-1)+ 
      labs(x="NMDS1",y="NMDS2",title=title) +
      annotate("text",label=stress,x=max(nmds.points$MDS1)*.9,y=max(nmds.points$MDS2),size=5)+
      geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))

    nmds.plot.nonmaes<-ggplot(nmds.points,aes(MDS1,MDS2)) + 
      geom_point(size=3,aes(shape=groups,color=groups)) + 
      labs(x="NMDS1",y="NMDS2",title=title) +
      annotate("text",label=stress,x=max(nmds.points$MDS1)*.9,y=max(nmds.points$MDS2),size=5)+
      geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
 } 

  ggsave(paste(outname,"NMDS.png",sep="."),nmds.plot,width=10,height=8)
  ggsave(paste(outname,"NMDS.pdf",sep="."),nmds.plot,width=10,height=8)
 
  ggsave(paste(outname,"NMDS.nonmaes.png",sep="."),nmds.plot.nonmaes,width=10,height=8)
  ggsave(paste(outname,"NMDS.nonames.pdf",sep="."),nmds.plot.nonmaes,width=10,height=8)  
}

nmds_plot(wu,"weighted_unifrac",title="Weighted unifrac NMDS")
nmds_plot(uw,"unweighted_unifrac",title="Unweighted unifrac NMDS")
