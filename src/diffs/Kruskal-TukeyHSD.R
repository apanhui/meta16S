args <- commandArgs(T)
library(vegan)
library(ggplot2)
File <- read.table(args[1],header=F)
Dirname <- args[2]

KruskalData <- matrix(NA,ncol=2,nrow=nrow(File))
colnames(KruskalData) <- c("Group","Pvalue")

for (i in 1:nrow(File)){
	file <- as.character(File[i,1])

	Rwdata <- read.table(file,header=T,check.names=F)

	Kruskal <- kruskal.test(Index~Group,data=Rwdata)
	AOV <- aov(Index~Group,data=Rwdata)
	HSD <- TukeyHSD(AOV)
	
	fileName <- basename(file)
	DiffName <- unlist( strsplit(fileName,split=".",fixed=T) )[1]
	KruskalData[i,] <- c(DiffName,Kruskal$p.value)
	HSD <- HSD[[1]]
	colnames <- c("Group",colnames(HSD))
	HSD2 <- cbind(rownames(HSD),HSD)
	colnames(HSD2)<- colnames
	write.table(HSD2,file=paste(Dirname,"/",DiffName,".TukeyHSD.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
#### draw boxplot
	GG <- ggplot(Rwdata,aes(Group,Index,color=Group)) +
	    geom_boxplot() +
	    geom_jitter(width=.2) +
	    labs(title=paste("the boxplot of ",DiffName,sep="")) +
	    theme(plot.title = element_text(size = 16)) +
	    guides(color=FALSE)

	ggsave(paste(Dirname,"/",DiffName,".boxplot.png",sep=""))
	ggsave(paste(Dirname,"/",DiffName,".boxplot.pdf",sep=""))
}

write.table(KruskalData,file=paste(Dirname,"/Kruskal.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
