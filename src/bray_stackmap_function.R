library(vegan)
library(ape)
args<-commandArgs(T)
data=read.table(args[1],sep="\t",head=T,row.names=1)
if(nrow(data)<2)
{
	cat("\n###\nthe taxa or otu count less than 2, can draw bray,but ugly , so quit\n###\n")
	quit("no")
}

dis <- vegdist(t(data),method="bray")
# dis
clus <- hclust(dis,"complete")
# clus$height
clus$height=clus$height*100
# ncol(data)

col=c("red","green3","blue","cyan","magenta","yellow","gray","purple","brown","darkred","darkmagenta","mediumvioletred","peru","lightcoral","darksalmon","sandybrown","black","darkgreen","forestgreen","mediumseagreen","limegreen","darkkhaki","mediumspringgreen","lawngreen","lightgreen","greenyellow","navy","midnightblue","darkviolet","darkslateblue","blueviolet","dimgray","lightslateblue","lightslategray","cornflowerblue","lightseagreen","violet","darkgray","darkturquoise","mediumturquoise","lightgray","powderblue","gainsboro","palegoldenrod","moccasin","lavender","blanchedalmond","papayawhip","linen","beige","oldlace","aliceblue","lightgoldenrodyellow","ghostwhite","floralwhite","mintcream")


a1=read.table(args[2],head=T,row.names=1,sep="\t")
stack_data = matrix(0,nrow = nrow(a1), ncol = ncol(a1))
for(i in 1:ncol(a1))
{
	stack_data[,i] = a1[,clus$order[i]]/sum(a1[,clus$order[i]])
}
clus$labels
clus$order
samples=rep("a",ncol(a1))
for(i in 1:length(samples))
{
	samples[i]=clus$labels[clus$order[i]]
}
# samples
#colnames(stack_data) <- levels(reorder(clus$labels,clus$order))
colnames(stack_data) <- samples
rownames(stack_data) <- rownames(a1)
#stack_data
b = rownames(stack_data)

### plot bray
#png(file="2.png")
#plot(clus,hang=-1,ylab = "Simlarity",axes=F)
#aa=axis(2,at=seq(0,100,by = 20),labels=rev(seq(0,100,by = 20)),las=2)
#dev.off()

ddd=as.phylo(clus)
edge_color=rep(1,nrow(ddd$edge))
sample_color<-rep(1,ncol(data))
mark=1     ## dont draw branch color
if(length(args)>=4)
{
	colms<-strsplit(args[4],",")
	classify<-as.numeric(colms[[1]])
	k=1
	j=1
	if(mark!=1)
	{
		for(i in 1:nrow(ddd$edge))
		{
			if(ddd$edge[i,2]!=classify[j])
			{
				edge_color[i]=col[k]
			}else{
				edge_color[i]=col[k]
				k=k+1
				j=j+1
			}
		}
	}
	sample_count<-length(ddd$tip.label)
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

mat=t(matrix(c(1,2),nrow=2))
node.depth.edgelength(ddd)
xrange=node.depth.edgelength(ddd)[1]

pdf(paste(args[3],".bray_stack.pdf",sep=""),16,8)
layout(mat,widths=c(2,3))
par(mar=c(7,3,3,0.8))
plot.phylo(ddd,label.offset=1,edge.width=2,edge.color=edge_color,tip.color=sample_color,y.lim=c(1-0.45,ncol(data)+0.45),show.tip.label=F)
text(x=max(xrange)*1.05,y=1:ncol(data),labels=samples,xpd=T)
axis(1,at=seq((xrange-50),xrange,by=10),labels=seq(0,100,by = 20))

print("------------------")
seq((xrange-50),xrange,by=10)
seq(0,100,by = 20)
print("------------------")

text(x=xrange/2,y=0.55-(ncol(data))/7,labels="Simlarity",xpd=T)
par(mar=c(7,0.8,3,as.numeric(args[4])))
dd=barplot(as.matrix(stack_data),xlab="Relative Abundance",main="",xlim=c(0,1),col=col[1:nrow(stack_data)],las=1,horiz=T,axisnames=F)
legend(x=1.03,y=median(dd),yjust=0.5,legend=b,cex=1.1,pt.cex=2.5,col=col[1:nrow(stack_data)],pch=15,xpd=T)
#legend(x=1.03,y=max(dd)+dd[1],legend=rev(b),cex=1.1,pt.cex=2.5,col=rev(col[1:nrow(stack_data)]),pch=15,xpd=T)
dev.off()


png(paste(args[3],".bray_stack.png",sep=""),width=1600,height=800)
layout(mat,widths=c(2,3,lcm(as.numeric(args[4]))))
par(mar=c(7,3,3,0.8))
plot.phylo(ddd,label.offset=1,edge.width=2,edge.color=edge_color,tip.color=sample_color,y.lim=c(1-0.45,ncol(data)+0.45),show.tip.label=F)
text(x=max(xrange)*1.05,y=1:ncol(data),labels=samples,xpd=T)
axis(1,at=seq((xrange-50),xrange,by = 10),labels=seq(0,100,by = 20))
text(x=xrange/2,y=0.55-(ncol(data))/7,labels="Simlarity",xpd=T)
par(mar=c(7,0.8,3,as.numeric(args[4])))
dd=barplot(as.matrix(stack_data),xlab="Relative Abundance",main="",xlim=c(0,1),col=col[1:nrow(stack_data)],las=1,horiz=T,axisnames=F)
dd
legend(x=1.03,y=median(dd),yjust=0.5,legend=b,cex=1.1,pt.cex=2.5,col=col[1:nrow(stack_data)],pch=15,xpd=T)
dev.off()

