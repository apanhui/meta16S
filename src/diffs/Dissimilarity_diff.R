#!/Bio/bin/Rscript

# Contact : ffhuang@genedenovo.com
args<-commandArgs(T)
if (length(args) != 3){
	usage = "
Version: 1.0

Descs:   significant difference between two or more groups of sampling units(adonis、anosim)

Usage:   Rscript Dissimilarity_diff.R <unweighted_unifrac.tsv,weighted_unifrac.tsv> <all.groups> <compare>

\n"

cat(usage)
quit("no")
}

# import packages
library(ggplot2)
library(vegan)
library(stringr)

# read data
dist <- read.table(args[1], header=T, row.names=1, check.names=F, "\t")
group <- read.table(args[2], sep="\t")

# get group information
split <- function(x){ y <- unlist( strsplit(x, "=") ) }
tmp <- apply(group, 1, split)

# remove blank
tmp <- str_trim(tmp) 
tmp <- matrix(tmp,ncol=nrow(group))
tmp <- t(tmp)

split <- function(x){ y <- unlist( strsplit(x, ",") ) }
Group_List <- tapply(tmp[,2], tmp[,1], split)
#Group_List

# 
diff <- args[3]
pairwises = unlist(strsplit(diff,","))
len <- length(pairwises)*2
results <- data.frame(diffs<-rep(NA,length=len), method<-rep(NA,length=len), pvalue<-rep(NA,length=len))
for (i in pairwises) {

	compares = unlist(strsplit(i,"&"))
	compares <- as.vector(compares)
	
	# get data
	ex <- function(x,Group_List){y<-Group_List[[x]]}
	part <- unlist(lapply(compares, ex, Group_List))
	part <- as.vector(part)
	p <- dist[part, part]
	p.dist <- as.dist(p)

	# get group information
	tmp <- function(x,Group_List){len <- length( Group_List[[x]]);y <-rep(x,len) }
	Group <- unlist(lapply(compares, tmp, Group_List))

	Ansoim <- anosim(p.dist,Group)

	Adonis <- adonis(p.dist~Group)

	index <- which(pairwises==i)
	i <- gsub("&","-vs-",i)
	results[2*index-1,] <- c(i,"anosim",Ansoim$signif)
	results[2*index,] <- c(i,"adonis",Adonis$aov.tab[1,6])
}

args[1] <- basename(args[1])
colnames(results)<-c("diffs","method","pvalue")
write.table(results,paste(basename(args[1]),".ssim.diff.xls",sep=""),row.names=F,sep="\t",quote=F)
