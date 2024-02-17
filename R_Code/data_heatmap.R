library(pheatmap)
library(tidyverse)
library(ggplotify)
library(heatmaply)

args = commandArgs(trailingOnly=TRUE)

width_ratio=10
height_ratio=8
margin=200

traits<-read.csv(args[1],header=TRUE,row.names=1, check.names=FALSE)
png(args[2],width=width_ratio*length(traits)+margin,height=height_ratio*length(traits[,1])+margin)
pheatmap(traits,  main=args[3])
dev.off()

