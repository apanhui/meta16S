#!/Bio/bin/Rscript
# name:    reads_stack.R
# func:    draw reads deal stack bar diagram for all samples 
# date:    2017-04-11
# author:  aipeng (pai@genedenovo.com)
# version: 1.0

library(ggplot2)

args<-commandArgs(T)

if (length(args) < 1)
{
	usage = "
Usage: Rscript reads_stack.R <all.data_process_stat.forFig.xls> <outdir>
"
	cat(usage)
	quit("no")
}

data = read.table(args[1],sep="\t",header=TRUE)

# calc the save fig size 
sample_num = length( unique(data$Sample) )
width  = 8
height = 6
if (sample_num <= 10) {
	width = 6
} else if (sample_num <= 30) {
	width = 8
} else if (sample_num >= 100) { 
	width = 16
} else {
	width = 8 + 8*(sample_num-30)/70 
}

# creat the color for 3D plot.
# the function to create default ggplot2 colors
gg_color_hue <- function(n) {                   
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# reverse the fill order 
level_factor  = factor(data$Reads,levels=unique(data$Reads),order=TRUE)
sample_factor = factor(data$Sample,levels=unique(data$Sample),order=TRUE)

num = length(unique(data$Reads))
#mycolors = rev(gg_color_hue(num))
mycolors = rev(c("#2ECDB3","#FEF104","#F58225","#A3218E","#EE1C27"))

ggplot(data,aes(sample_factor,Number)) + 
geom_bar(aes(fill=level_factor),stat="identity") + 
labs(x="Sample Names",y="Number of Reads/Tags",fill="") + 
scale_fill_manual(values=mycolors) + 
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

ggsave(paste(args[2],"reads_stack.pdf",sep="/"),width=width,height=height)
ggsave(paste(args[2],"reads_stack.png",sep="/"),width=width,height=height)

ggplot(data,aes(sample_factor,Number)) + 
geom_bar(aes(fill=level_factor),stat="identity",position="fill") + 
labs(x="Sample Names",y="Percentage of All Reads",fill="") + 
scale_fill_manual(values=mycolors) + 
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

ggsave(paste(args[2],"reads_fill.pdf",sep="/"),width=width,height=height)
ggsave(paste(args[2],"reads_fill.png",sep="/"),width=width,height=height)
