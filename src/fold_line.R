#!/Bio/bin/Rscript
# name    : fold_line.R
# func    : draw fold lines for alpha diversity
# version : 1.0
# date    : 2017-06-05
# contact : pai@genedenovo.com
#-------------------------------------------------
args<-commandArgs(T) 

if (length(args) != 1){
  usage = "
Version: 1.0

Descs  : An Script for draw alpha diversity fold lines

Usage  : Rscript fold_line.R <values.txt>
\n"
 cat(usage)
 quit("no")
}

# import packages
library(ggplot2)

outname = unlist(strsplit(args[1],"./"))
outname = gsub(pattern=".forGGPLOT2.txt",replacement="",outname[length(outname)])

ylab = ""
if (length(grep("observed_species",outname)) == 1){
  ylab = "Number of OTUs observed"
} else if (length(grep("shannon",outname)) == 1) {
  ylab = "Shannon Index"
} else if (length(grep("chao1",outname)) == 1) {
  ylab = "Chao1 Index"
} else {
  ylab = outname
}

data = read.table(args[1],header=TRUE,check.names=FALSE,sep="\t",comment.char="")

# read the orignal order of samples/groups
sample_factor = factor(data[[1]],levels=unique(data[[1]]),order=TRUE)

# read the legend title, the first filed name (Sample/Group)
colnames = colnames(data)
legend_title = colnames[1]

# read the colors 
# the function to create default ggplot2 colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

mycolors = levels(data$color)
num = length(mycolors)
sample_num = length(unique(data[[1]]))
if (num != sample_num)
{
  mycolors = gg_color_hue(sample_num)
}

# init the legend.position
mytheme = theme()
if (num <= 10 && length(grep("observed_species",outname)) == 1){
  mytheme = theme(legend.position = c(0,1), legend.justification = c("left", "top"))
} else if (num <= 25) {
  mytheme = theme(legend.position = "right")
} else {
  mytheme = theme(legend.position = "none")
}

data$ymin = data$val - data$se
data$ymax = data$val + data$se

ggplot(data) + 
  geom_line(aes(x=x,y=val,color=sample_factor),size=1.2) + 
  #geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax,group=sample_factor,color=sample_factor),size=1) + 
  scale_color_manual(values = mycolors) + 
  labs(x="Number of Tags Sampled",y=ylab,color=legend_title) + 
  theme_bw() + mytheme

ggsave(paste(outname[1],"png",sep="."))
ggsave(paste(outname[1],"pdf",sep="."))
