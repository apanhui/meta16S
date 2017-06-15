args <- commandArgs(T)
library(vegan)
library(ggplot2)
File <- read.table(args[1],header=F)
Dirname <- args[2]

TtestData <- matrix(NA,ncol=2,nrow=nrow(File))
colnames(TtestData) <- c("Group","Pvalue")
WilcoxData <- matrix(NA,ncol=2,nrow=nrow(File))
colnames(WilcoxData) <- c("Group","Pvalue")

for (i in 1:nrow(File)){
	file <- as.character(File[i,1])

	Rwdata <- read.table(file,header=T,check.names=F)

	Ttest <- t.test(Index~Group,data=Rwdata)
	Wilcox <-wilcox.test(Index~Group,data=Rwdata)

	fileName <- basename(file)
	DiffName <- unlist( strsplit(fileName,split=".",fixed=T) )[1]
	TtestData[i,] <- c(DiffName,Ttest$p.value)
	WilcoxData[i,] <- c(DiffName,Wilcox$p.value)
##### draw boxplot
#	png(file=paste(Dirname,"/",DiffName,".boxplot.png",sep=""),width=800,height=600)
#	pdf(file=paste(Dirname,"/",DiffName,".boxplot.pdf",sep=""),width=8,height=6)
	GG <- ggplot(Rwdata,aes(Group,Index,color=Group)) +
	    geom_boxplot() +
	    geom_jitter(width=.2) +
	    labs(title=paste("the boxplot of ",DiffName,sep="")) +
	    theme(plot.title = element_text(size = 16)) +
	    guides(color=FALSE)

	ggsave(paste(Dirname,"/",DiffName,".boxplot.png",sep=""))
	ggsave(paste(Dirname,"/",DiffName,".boxplot.pdf",sep=""))
}

write.table(TtestData,file=paste(Dirname,"/Ttest.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
write.table(WilcoxData,file=paste(Dirname,"/Wilcox.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)



