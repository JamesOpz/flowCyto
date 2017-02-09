# JO Initial rough script 30/01/17
# load the SPADE library
library(spade)
library(gplots)
library(pheatmap)
library(RColorBrewer)

# To Do:
# In theory we want to split the populations by keyword and add an pseudovalue column so that we can map the cellular idenitites back on the nodes

# Use a file exported from FlowJo, must contain the t-SNE values within it
# Should export compensated paramters
data_file_path = '942_UNCOMP_TAM Ter119 perpc55 Cd45 FITC F480 BV410 MHC II PE CD11c 780 CD11b 510_942 cd206 A_CD45+.fcs'

# QC
fsc_file <- flowCore::read.FCS(data_file_path)
fsc_file

# Identify number if cell clusters, not 100% sure how many is optimal, but it appears that overclustring is important
# Clustering should be proportional to the cell number as you dont want loads of nodes with no cells in it 
Clusters <- 100

# Name of output directory  
output_dir <- file.path("942_TAM_ALLPARAMS")

# Select the markers to run the clustering on, this is going to be the tSNE channels
markers <- c('APC-A', 'PE-A', 'Pacific Blue-A', 'AmCyan-A', 'APC-Cy7-A')

# not sure how to tramsform the other data, i don't think it really matters
transforms <- c('APC-A'=flowCore::arcsinhTransform(a=0, b=0.01), 'AmCyan-A'=flowCore::arcsinhTransform(a=0, b=0.01), 'PE-A'=flowCore::arcsinhTransform(a=0, b=0.01), 'Pacific Blue-A'=flowCore::arcsinhTransform(a=0, b=0.01), 'APC-Cy7-A'=flowCore::arcsinhTransform(a=0, b=0.01))

#Execute core SPADE algorithm
SPADE.driver(data_file_path, out_dir=output_dir, transforms=transforms, cluster_cols=markers,k= Clusters)

# Look at the files generated
grep("^libloc",dir(output_dir),invert=TRUE,value=TRUE)

# loade the graph object into memory
mst_graph <- igraph:::read.graph(paste(output_dir,"mst.gml",sep=.Platform$file.sep),format="gml")

# Load the node marker intesity values for heatmapping 
layout <- read.table(paste(output_dir,"layout.table",sep=.Platform$file.sep))

# Plot all the trees!! (In kamada kawai force driected graph mapping) 
SPADE.plot.trees(mst_graph, output_dir, out_dir=paste(output_dir,"pdf",sep=.Platform$file.sep), layout=as.matrix(layout))

# Heatmap the sample table

# Import the node median fluoresencs values
table_file_path <- paste(data_file_path, '.density.fcs.cluster.fcs.anno.Rsave_table.csv', sep='')
table_path <- file.path(output_dir,'tables', 'bySample', table_file_path)
table_path <- normalizePath(table_path)
SPADE_nodes <- as.matrix(read.csv(table_path))

# Prune the table values so just the ones we want
# Can automate this
# removed count for now
cluster_panel <- c('ID', 'mediansAPC.A_clust', 'mediansPE.A_clust', 'mediansPacific.Blue.A_clust', 'mediansAmCyan.A_clust', 'mediansAPC.Cy7.A_clust')
prunedByChannel_SPADE_nodes <- SPADE_nodes[,cluster_panel]

# prune out minimum channels, do this later
# pruned_SPADE_nodes <- prunedByChannel_SPADE_nodes[which(prunedByChannel_SPADE_nodes[,'count'] > 5)]

# prune the table IDs into rownames
# This is stupid as its dependant on the number of markers used
x <- length(cluster_panel)
rownames(prunedByChannel_SPADE_nodes) <- prunedByChannel_SPADE_nodes[,1]
prunedByChannel_SPADE_nodes <- prunedByChannel_SPADE_nodes[,2:x]

# arcsinh width basis is fidly, 0.01 seems to work quite well
# If the data is previously transformed by SPADE no need
# arcsinh_NODES <- asinh(prunedByChannel_SPADE_nodes*0.01) 
arcsinh_NODES <- prunedByChannel_SPADE_nodes

#Rename the columns
colnames(arcsinh_NODES) = c('CD206', 'MHCII', 'F4/80', 'CD11b', 'CD11c')

# Do the heatmapping
#Default heat scale- black, yellow, blue-- may need to change range depending on data values
heat_palette <- colorRampPalette(c("black", "yellow", "#FAF7C9"))
hot_pal <- heat_palette(80)
bluecols <- brewer.pal(80, 'Blues')
blue_pallette <- colorRampPalette(bluecols)
ncols <- 81
bluecols2 <- newcol(ncols)

#heat_palette break for ramp
heat_pairs.breaks <- seq(-1, 7, by=0.1)

# Blues breaks - note for pheatmap breaks must be a seq covering the range of values in x 1 element longer than colour vector.
blue_break <- seq(0, 4,by=0.5)

# Need to fix messed up formatting of dendrogram

# dev.new(width=10, height=20)
pdf('pheatmap_ALLPARAMS.pdf', paper='a4', width = 12, height = 15)

# Heatmap.2 execution
# heatmap.2(arcsinh_NODES, main="SPADE node medians",breaks=scaled_pairs.breaks, dendrogram="both",trace="none",Rowv = TRUE, Colv=TRUE, revC=FALSE,symm=FALSE, symkey=FALSE,symbreaks=FALSE,scale="none",cexRow=0.6,cexCol=1,key=TRUE,col=heat_palette, density.info="none", margins = c(8,3), lmat = matrix(c(4,2,3,1), nrow=2, ncol=2) , lhei = c(0.1*2,0.9*2), lwid = c(0.3*(1/50),0.7*(1/50)))

# For pheatmap, calculate the optimal fontsize to fit the row
fontsize_row = 12 - nrow(arcsinh_NODES) / 15

# pheatmap execution
pheatmap(arcsinh_NODES, color = hot_pal, breaks = heat_pairs.breaks,  scale = 'none', cellwidth = 9, cellheight =  7, fontsize_row = fontsize_row)

dev.off()
