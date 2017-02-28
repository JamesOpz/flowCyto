# JO working script with updates from Shahram
# 15/02/2017

# To Do list:
# Can I split the populations by keyword 
# Then add an pseudovalue column so that we can map the cellular idenitites back on the nodes



# load the SPADE library
library(spade)

# load flowR functional libraries
library(openCyto)
library(ggcyto)

# Load the heatmap functionality libraries
# Note heat map will transition to Superheat  
library(gplots)
library(pheatmap)
library(superheat)
library(RColorBrewer)


# Use a file exported from FlowJo, must contain the t-SNE values within it
# Should export uncompensated parameters


data_file_path = '1000_TAM Ter119 perpc55 Cd45 780 F480 BV410 MHC II PE Ly6C 488 CD11b 510_1007 CD206A_Single Cells.fcs'

QC_fcs_file <- function(data_file_path) {
  # QC with user interface to check
  fsc_file <- flowCore::read.FCS(data_file_path)
  fsc_file

}

# This function should make dotplots of the tSNE graphs
# tSNE graphs should be heat cloloured by channel
# With scale transform specified

plot_tSNE <- function(data_file_path, channel=NULL, transform=NULL) {}


# Identify number if cell clusters, not 100% sure how many is optimal, 
# Detect number of clusters with NBclust?
# it appears that overclustring is important
# Clustering should be proportional to the cell number as you dont want loads of nodes with no cells in it 

RUN_SPADE_tSNE <- function(data_file_path, clusters=100, n=3, transforms='flowCore::arcsinhTransform(a=0, b=0.01))') {
  Clusters <- 100

  # Name of output directory
  out <- strsplit(data_file_path,' ')
  out_path <- out[[1]][1]
  output_dir <- file.path(out_path)

  # Creat master directory to put repeats in
  
  # Select SPADE to run on the tSNE channels the
  # Note that the tSNE settings are embedded into channel name and must be the same  
  markers <- c('tSNE_X_P_30_E_200_I_1000_T_0.2', 'tSNE_Y_P_30_E_200_I_1000_T_0.2')

  # This is for specific marker transforms
  # transforms <- c('APC-A'=flowCore::arcsinhTransform(a=0, b=0.01), 'AmCyan-A'=flowCore::arcsinhTransform(a=0, b=0.01), 'PE-A'=flowCore::arcsinhTransform(a=0, b=0.01), 'Pacific Blue-A'=flowCore::arcsinhTransform(a=0, b=0.01), 'APC-Cy7-A'=flowCore::arcsinhTransform(a=0, b=0.01))

  # Execute core SPADE algorithm in replicates of 3
  # To generate 3 cluster tables to average and then heatmap
  
  for(i in 1:n){
    
    output_dir_n <- paste(output_dir, i, sep = '_') 
    
    # execute core SPADE command
    spade::SPADE.driver(data_file_path, out_dir=output_dir_n, transforms=transforms, cluster_cols=markers,k= Clusters)

    # load the graph object into memory
    mst_graph <- igraph:::read.graph(paste(output_dir_n,"mst.gml",sep=.Platform$file.sep),format="gml")

    # Load the node marker intesity values for heatmapping 
    layout <- read.table(paste(output_dir_n,"layout.table",sep=.Platform$file.sep))

    # Plot all the trees!! (In kamada kawai force driected graph mapping) 
    SPADE.plot.trees(mst_graph, output_dir_n, out_dir=paste(output_dir_n,"pdf",sep=.Platform$file.sep), layout=as.matrix(layout))

  }
}

# Include markers as vector to determine cluster markers
# This is currently set up on TAM panel
av_table_heatmap_SPADE <- function(data_file_path, n=3) {

  # Function to average the three tables
  # Then write out the average table and heatmap it

  out <- strsplit(data_file_path,' ')
  out_path <- out[[1]][1]
  output_dir <- file.path(out_path)  

  for(p in 1:n) {
  
    output_dir_n <- paste(output_dir, p, sep = '_')  

    # Import the node median fluoresencs values and add them together
    table_file_path <- paste(data_file_path, '.density.fcs.cluster.fcs.anno.Rsave_table.csv', sep='')
    table_path <- file.path(output_dir_n,'tables', 'bySample', table_file_path)
    table_path <- normalizePath(table_path)
    SPADE_nodes <- as.matrix(read.csv(table_path))
  
  if(p == 1) {
    
    combined_SPADE_nodes <- SPADE_nodes
    print(SPADE_nodes[5,9])
  }
  
  if(p > 1) {
  
    combined_SPADE_nodes <- combined_SPADE_nodes + SPADE_nodes
    print(SPADE_nodes[5,9])  
  }
  }

  #average the n combined tables  
  av_SPADE_nodes <- combined_SPADE_nodes / n
  print(n)
  print(combined_SPADE_nodes[5,9])
  print(av_SPADE_nodes[5,9])
  
  # Prune the table values so just the ones we want
  # Use markers from input
  
  cluster_panel <- c('ID', 'mediansX.APC.A', 'mediansX.PE.A', 'mediansX.Pacific.Blue.A', 'mediansX.FITC.A','mediansX.APC.Cy7.A', 'mediansX.AmCyan.A')
  prunedByChannel_SPADE_nodes <- av_SPADE_nodes[,cluster_panel]
  
  # write pruned by channel nodes to excel
  write.csv(prunedByChannel_SPADE_nodes, file = 'av_spade_nodes.csv', sep = ",")
  
  # prune out minimum clusters at 1%, do this later
  # pruned_SPADE_nodes <- prunedByChannel_SPADE_nodes[which(prunedByChannel_SPADE_nodes[,'count'] > 5)]

  # prune the table IDs into rownames
  # This is stupid as its dependant on the number of markers used
  x <- length(cluster_panel)
  rownames(prunedByChannel_SPADE_nodes) <- prunedByChannel_SPADE_nodes[,1]
  prunedByChannel_SPADE_nodes <- prunedByChannel_SPADE_nodes[,2:x]

  # arcsinh width basis is fidly, 0.01 seems to work quite well
  # If the data is previously transformed by SPADE no need
  arcsinh_NODES <- asinh(prunedByChannel_SPADE_nodes*0.007) 
  # arcsinh_NODES <- prunedByChannel_SPADE_nodes
  
  # Shift to row z-score ot bet represent the relative protein expression
  
  scale_nodes <- scale(prunedByChannel_SPADE_nodes) 

  #Rename the columns - this should be done with a double vector of markers and stains
  colnames(arcsinh_NODES) = c('CD206', 'MHCII', 'F4_80', 'Ly6C', 'CD45', 'CD11b')
  
  
  
  # print(prunedByChannel_SPADE_nodes)
  
  # print(arcsinh_NODES) 
  
  # print(scale_nodes)
  
  # Do the heatmapping
  
  #Default heat scale- black, yellow, blue-- may need to change range depending on data values
  heat_palette <- colorRampPalette(c("black", "yellow", "#FAF7C9"))
  hot_pal <- heat_palette(80)

  # heat_palette break for ramp
  # heat_pairs.breaks <- seq(-1, 7, by=0.1)

  # Red blue heatmap
  
  print(max(arcsinh_NODES))
  print(min(arcsinh_NODES))
  
  cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  blue_breaks <- seq(min(prunedByChannel_SPADE_nodes), max(prunedByChannel_SPADE_nodes), by = (max(prunedByChannel_SPADE_nodes) - min(prunedByChannel_SPADE_nodes)) / 256) 
  
  # Need to fix messed up formatting of dendrogram 
  # shifting to superheat package
  
  pdf_outfile <- paste(output_dir, '.pdf', sep = '')
  
  pdf(pdf_outfile, paper='a4', width = 12, height = 15)

  # For pheatmap, calculate the optimal fontsize to fit the row
  fontsize_row = 12 - nrow(arcsinh_NODES) / 15

  # pheatmap execution
  pheatmap(arcsinh_NODES, color = hot_pal,  scale = 'none', cellwidth = 9, cellheight =  7, fontsize_row = fontsize_row)

  # superheat execution
  
  dev.off()
  
}