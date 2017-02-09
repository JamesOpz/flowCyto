#Heatmaps for SPADE node medians

#Set working directory to folder where text files exported from SPADE are located
setwd(" ")

#Install package gplots of not already in library
#Load gplots
library(gplots)

#Import medians from tab-delimited text file (from SPADE)
SPADE_nodes = as.matrix(abs(read.table("file.txt",sep="\t",header=TRUE)))

#Default heat scale- black, yellow, blue-- may need to change range depending on data values
heat_palette <- colorRampPalette(c("black", "yellow", "#FAF7C9"))
pairs.breaks <- c(seq(0,1, by=0.01),seq(1.01, 2, by=0.01),seq(2.01,3, by=0.01),seq(3.01,4, by=0.01),seq(4.01, 5, by=0.01),seq(5.01,6, by=0.01))

#Generate heatmap and dendrogram
heatmap.2(SPADE_nodes, main="SPADE node medians",breaks=pairs.breaks,dendrogram="both",trace="none",Rowv = TRUE, Colv=TRUE, revC=FALSE,symm=FALSE, symkey=FALSE,symbreaks=FALSE,scale="none",cexRow=1,cexCol=0.8,key=TRUE,col=heat_palette, density.info="none")

