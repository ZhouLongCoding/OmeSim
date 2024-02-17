library(ggraph)
library(igraph)

args = commandArgs(trailingOnly=TRUE)


network=read.csv(args[1])
g <- graph_from_data_frame(network, directed = TRUE)

#corr_data<-read.csv(args[1],header=TRUE,row.names=1, check.names=FALSE)
png(args[2],width=args[3],height=args[3])
plot(g, edge.arrow.size = 0.5, vertex.label.color = "blue", vertex.size = 8)
title(main = args[4], col.main = "black", font.main = 20)
dev.off()
