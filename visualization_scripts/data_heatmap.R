library(pheatmap)
library(tidyverse)
library(ggplotify)
library(heatmaply)

args = commandArgs(trailingOnly=TRUE)


traits <- read.csv(args[1], header=TRUE, check.names=FALSE)
colnames(traits)[1] <- "Gene_ID"
traits <- aggregate( . ~ Gene_ID,data=traits, max)
rownames(traits) <- traits[,1]
traits <- traits[,-1]


if (ncol(traits) <= 1000) {
  width_ratio <- 10
  height_ratio <- 8
  margin <- 200
  
  png(args[2],
      width  = width_ratio * ncol(traits) + margin,
      height = height_ratio * nrow(traits) + margin)
  
  pheatmap(traits,  
           main = args[3],
           fontsize = 20,
           fontsize_row = 10,
           fontsize_col = 10)
  dev.off()
  
} else {
  width_ratio <- 5
  height_ratio <- 5
  margin <- 100
  
  png(args[2],
      width  = width_ratio * ncol(traits) + margin,
      height = height_ratio * nrow(traits) + margin)
  
  pheatmap(traits,  
           main = args[3],
           fontsize = 20,
           fontsize_row = 5,
           fontsize_col = 5)
  dev.off()
}


