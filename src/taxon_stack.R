#!/Bio/bin/Rscript
# name:    taxon_stack.R
# func:    draw taxonomy stack bar diagram for all samples 
# date:    2017-04-11
# author:  aipeng (pai@genedenovo.com)
# version: 1.0

library(ggplot2)

args<-commandArgs(T)

if (length(args) < 1)
{
	usage = "
Usage: Rscript taxon_stack.R <all.taxonomy.stat.forFig.xls>
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

# reverse the fill order 
level_factor  = factor(data$Level,levels=unique(data$Level),order=TRUE)
sample_factor = factor(data$Sample,levels=unique(data$Sample),order=TRUE)

ggplot(data,aes(sample_factor,Number)) + 
geom_bar(aes(fill=level_factor),stat="identity") + 
labs(x="Sample Names",y="Number of Tags",fill="Taxonomy") + 
#guides(fill = guide_legend(reverse=TRUE)) + 
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

ggsave("taxonomy_stack.pdf",width=width,height=height)
ggsave("taxonomy_stack.png",width=width,height=height)

ggplot(data,aes(sample_factor,Number)) + 
geom_bar(aes(fill=level_factor),stat="identity",position="fill") + 
labs(x="Sample Names",y="Percentage of All Tags",fill="Taxonomy") + 
#guides(fill = guide_legend(reverse=TRUE)) + 
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

ggsave("taxonomy_fill.pdf",width=width,height=height)
ggsave("taxonomy_fill.png",width=width,height=height)
