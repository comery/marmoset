#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

library(dplyr)
library(ggplot2)
library(reshape2)
library(qgraph)
library(tidyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

pairs<-(read.csv("pairs.noAutosomes.noX.inter.list", sep = "\t")[, -c(9,10)])
chr.sizes<-read.csv("chrs.txt", sep = "\t", header = FALSE)

ass<- pairs%>% select(chrom1, chrom2)

ass$chrom1 = factor(ass$chrom1,levels=sort(union(ass$chrom1,ass$chrom2)))
ass$chrom2 = factor(ass$chrom2,levels=sort(union(ass$chrom1,ass$chrom2)))
res = as.data.frame.matrix(table(ass))

counts<-data.frame(table(ass))

counts<-dplyr::left_join(counts, chr.sizes, by = c("chrom1" = "V1"))
counts<-dplyr::left_join(counts, chr.sizes, by = c("chrom2" = "V1"))

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
counts <- transform(counts, min = pmin(as.numeric.factor(V2.x), as.numeric.factor(V2.y)))

counts<-counts %>% mutate(ncount = Freq /min)
matrix<-counts %>% select(chrom1, chrom2, ncount)

matrix<-as.data.frame.matrix(dcast(matrix, chrom1~chrom2))
matrix <- data.frame(matrix[,-1], row.names=matrix[,1])

matrix[matrix<=0.001]<-0
matrix[matrix>0.0001]<-1


jpeg('distance graph1.jpg', width=10000, height=10000, unit='px')
g1<-graph_from_adjacency_matrix( as.matrix(matrix))
plot(g1,vertex.size=4, edge.arrow.size=0.5)
dev.off()

jpeg('distance graph2.jpg', width=10000, height=10000, unit='px')
qgraph(matrix, layout='circle', threshold=0.005, vsize=1, esize=1, arrows=FALSE, borders=FALSE)
dev.off()

jpeg('heatmap.jpg', width=10000, height=10000, unit='px')
heatmap(as.matrix(res))
dev.off()

