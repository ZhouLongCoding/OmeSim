library(pheatmap)
library(tidyverse)
library(ggplotify)
library(heatmaply)

args = commandArgs(trailingOnly=TRUE)

width_ratio=10
height_ratio=8
margin=200

corr_data <- read.csv(args[1], header=TRUE, check.names=FALSE)

if(length(unique(corr_data[,1]))==dim(corr_data)[1]){
  rownames(corr_data) <- corr_data[,1]
  corr_data <- corr_data[,-1]
}else{
  duplicated_columns <- which(duplicated(corr_data[,1]))
  corr_data <- corr_data[-duplicated_columns,-(duplicated_columns+1)]
  rownames(corr_data) <- corr_data[,1]
  corr_data <- corr_data[,-1]
}

png(args[2],width=width_ratio*length(corr_data)+margin,height=height_ratio*length(corr_data[,1])+margin)
pheatmap(corr_data, main=args[3], cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()
