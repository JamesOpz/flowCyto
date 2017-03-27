# JO working script with updates from Shahram
# 01/03/2017
# V2 - integrated functional script

# To Do list:
# Calculate CVs across nodes 
# Make the heatmap all fancy 
# Make the colour palette really nice
# Repaet channel poltting
# Re-base 
# Driver funtion

# Need to trace modify the write.FCS fucntion to alter the $PnR - done

viSNE.analysis.driver <- function(FCS_file_path, csv_file_path = NULL, covertFCS = FALSE, markers, fpanel){
  
  # Plot tSNE values for all channels - incomplete
  plot

  # Generate new FCS file if specified
  if(convertFCS = TRUE){
    
    read_csvToFcs(csv_file_path)
   
  }
  
  
  # QC the fcs file to log - incomplete
  
  # Determine number of clusters - incomplete
  
  
  # Run SPADE clusteirng 
  RUN_SPADE_tSNE(fcs_file)
  
  # Plot a pretty cluster_map
  plot_pretty_clusters(fcs_file)
  
  # Heatmap
  
  
  
}

# must find cluster file name
plot_pretty_clusters <- function(fcs_file_cluster){
  
  # Name of output directory
  out <- strsplit(fcs_file_cluster,' ')
  out_path <- out[[1]][1]
  out_path_1 <- paste(out_path, "1", sep = "_")
  output_dir <- file.path(out_path)
  
  # Find the cluster file in the file
  cluster_file <- fcs_file_cluster
  
  # read and display the cluater_fcs file
  cluster_fcs <- flowCore::read.FCS(cluster_file,truncate_max_range = FALSE)
  
  # extrat the data from within
  tSNE_clust <- exprs(cluster_fcs)
  
  clust_df <- data.frame(tSNE_clust)
  
  # define color schema
  # must be same colour length as the number clusters used for SPADE
  # Not this is in a 12 colour rainbow, must figure out how to get mroe colours
  clust_df$Col <- palette(rainbow(200))[as.numeric(cut(clust_df$cluster, breaks = 200))]
  
  #plot as before
  plot.default(clust_df$tSNE_X_P_30_E_200_I_1000_T_0.2, 
               clust_df$tSNE_Y_P_30_E_200_I_1000_T_0.2, 
               pch=19, 
               cex=0.3,
               col = clust_df$Col,
               xlab = 'tSNE 1',
               ylab = 'tSNE 2'
  )
  
}

# Find the optimum number of clusters in the data structure
# Currently non-functional
test_num_clusters <- function(fcs_file){
  
  # read and display the cluater_fcs file
  cluster_fcs <- flowCore::read.FCS(fcs_file,truncate_max_range = FALSE)
  
  # extrat the data from within
  tSNE_clust <- exprs(cluster_fcs)
  
  clust_df <- data.frame(tSNE_clust)
  
  tSNE_vals <- clust_df[, c('tSNE_X_P_30_E_200_I_1000_T_0.2', 'tSNE_Y_P_30_E_200_I_1000_T_0.2')]
  
  tSNE_vals_clust <- NbClust(tSNE_vals, min.nc = 50, max.nc = 300, method = 'complete') 
  
}

# Correctly annotate the heatmap with user defined annotations
# Generate new heatmap
# Currently does'nt work
annotate.heatmap <- function(){
  
  # Take the heatmap clusters and split them adding annotations
  
  
  
  
  if(is.list(annotations)){
    
    # generate vector for categorizing 
    clust_df$categories <- 'Un'
    for(l in 1:length(annotations)){
      
      positions <- which(clust_df$cluster %in% annotations[[l]])
      clust_df$categories[positions] <- names[l]
    }
    
    clust_df$categories <- as.factor(clust_df$categories)
    
    #create the new color pallet
    num_colors <- nlevels(clust_df$categories)
    coloring <- rainbow(num_colors)
    
    plot.default(clust_df$tSNE_X_P_30_E_200_I_1000_T_0.2, 
                 clust_df$tSNE_Y_P_30_E_200_I_1000_T_0.2, 
                 pch=19, 
                 cex=0.3,
                 col = coloring[clust_df$categories],
                 xlab = 'tSNE 1',
                 ylab = 'tSNE 2'
    )
    
    legend(
      x ="topleft",
      legend = paste(levels(clust_df$categories)), # for readability of legend
      col = coloring,
      pch = 19, # same as pch=20, just smaller
      cex = 0.7, # scale the legend to look attractively sized
      bty = 'n'
    )
    
  }
}

read_csvToFcs <- function(dat) {
  
  
  dta <- read.csv(dat)
  
  head(dta)
  
  # you need to prepare some metadata
  meta <- data.frame(name=dimnames(dta)[[2]],
                     desc=paste('this is column',dimnames(dta)[[2]],'from your CSV')
  )
  meta$range <- apply(apply(dta,2,range),2,diff)
  meta$minRange <- apply(dta,2,min)
  meta$maxRange <- apply(dta,2,max)
  
  head(meta)
  
  mat.dta <- as.matrix(dta)
  
  # all these are required for the following steps to work
  
  # a flowFrame is the internal representation of a FCS file
  ff <- new("flowFrame",
            exprs=mat.dta,
            parameters=AnnotatedDataFrame(meta)
  )
  
  # a simple plot to check that it worked
  dev.new(width=1.75, height=5)
  
  xyplot(SSC.A~FSC.A,ff)
  
  out <- strsplit(dat,'.csv')
  output_dir <- file.path(out)
  
  # now you can save it back to the filesystem
  # Need to find a permanent way to edit the filesystem
  write.FCS(ff,paste(output_dir,'FCS',sep='.'))
  
  
  
  
}

# load the ressescary libraries
load_spade_libs <- function() {
  
  # load the SPADE library
  library(spade)
  library(Biobase)
  library(flowViz)
  # load flowR functional libraries
  library(openCyto)
  library(ggcyto)
  
  # Load the heatmap functionality libraries
  # Note heat map will transition to Superheat  
  library(gplots)
  library(pheatmap)
  library(superheat)
  library(RCofillorBrewer)
  
}



QC_fcs_file <- function(data_file_path) {
  
  
  # QC with user interface to check
  fsc_file <- flowCore::read.FCS(data_file_path, truncate_max_range = FALSE)
  fsc_file

}

# This function should make dotplots of the tSNE graphs
# tSNE graphs should be heat cloloured by channel
# With scale transform specified use asinh with t factor = 0.0095

color.bar <- function(lut, pdf_file , min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  pdf(pdf_file)
  # dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }

  dev.off()
  
  }

plot_tSNE <- function(csv_file_path, transform=0.02) {
  
  # read the accense data in put and plot the flowJO tSNE points
  tSNE_data <- read.csv(csv_file_path, sep = ',')
  
  # Create the color ramp
  # spectralPal <- colorRampPalette(brewer.pal(11,"Spectral"))(100)
  rbPal <- colorRampPalette(c('blue', 'grey', 'red'))
  heatPal <- colorRamps::matlab.like(20)
  
  
  
  # Create heat colour scale and add it as a vector column
  for (k in 3:(ncol(tSNE_data)-2)) {
  
    tSNE <- as.matrix(tSNE_data[k])
    
    asinh_tSNE <- asinh(tSNE*transform)
    
    # Assign data to heat values
    tSNE_cols <- as.numeric(cut(asinh_tSNE, breaks = c(round(min(asinh_tSNE), 2), seq( -0.1, round(quantile(asinh_tSNE, 0.9995), 2), len = 99))))    
    plot_Cols <- rbPal(99)[tSNE_cols]
  
    # Create the output_file 
    pdf_outfile <- paste(colnames(tSNE_data[k]), 'asinh', transform , '.pdf', sep = '_')
    pdf(pdf_outfile, paper='a4', width = 50, height = 40)
    
    # split the pdf into two bits for legend
    # layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
    
    
    
    plot.default(tSNE_data$tSNE_X_P_30_E_200_I_1000_T_0.2, 
                tSNE_data$tSNE_Y_P_30_E_200_I_1000_T_0.2, 
                pch=19, 
                cex=0.45,
                col = plot_Cols,
                xlab = 'tSNE 1',
                ylab = 'tSNE 2'
                )
    
    dev.off()
    
    pdf_bar <- paste(colnames(tSNE_data[k]), 'Legend', '.pdf', sep = '_')
    
    color.bar(rbPal(80), pdf_bar, min = -0.1, max = round(quantile(asinh_tSNE, 0.9995), 2))
  
}

}

# Identify number if cell clusters, not 100% sure how many is optimal, 
# Detect number of clusters with NBclust?
# it appears that overclustring is important
# Clustering should be proportional to the cell number as you dont want loads of nodes with no cells in it 
# Does SPADE have to be run with the same random seed everytime for conistency?

RUN_SPADE_tSNE <- function(data_file_path, clusters=100, n=3, transforms='none', markers = c('tSNE_X_P_30_E_200_I_1000_T_0.2', 'tSNE_Y_P_30_E_200_I_1000_T_0.2'))  {
  

  # Name of output directory
  out <- strsplit(data_file_path,' ')
  out_path <- out[[1]][1]
  output_dir <- file.path(out_path)

  # Execute core SPADE algorithm in replicates of 3
  # To generate 3 cluster tables to average and then heatmap
  
  for(i in 1:n){
    
    output_dir_n <- paste(output_dir, i, sep = '_') 
    
    # execute core SPADE command
    # reset the seed for consistent values
    set.seed(1)
    spade::SPADE.driver(data_file_path, out_dir=output_dir_n, transforms='none', cluster_cols=markers, k=clusters)

    # load the graph object into memory
    mst_graph <- igraph:::read.graph(paste(output_dir_n,"mst.gml",sep=.Platform$file.sep),format="gml")

    # Load the node marker intesity values for heatmapping 
    layout <- read.table(paste(output_dir_n,"layout.table",sep=.Platform$file.sep))

    # Plot all the trees!! (In kamada kawai force driected graph mapping) 
    SPADE.plot.trees(mst_graph, output_dir_n, out_dir=paste(output_dir_n,"pdf",sep=.Platform$file.sep), layout=as.matrix(layout))

  }
}


# Include markers as vector to determine cluster markers
# This is currently set up on MSC panel

av_table_heatmap_SPADE <- function(data_file_path, n=3, transform = 0.02, annotations = NULL, clust_df = NULL, names = NULL, panel) {

  # Function to average the three tables
  # Then write out the average table and heatmap it

  out <- strsplit(data_file_path,' ')
  out_path <- out[[1]][1]
  output_dir <- file.path(out_path)  
  replicates <- list()
  count <- 1
  
  for(p in 1:n) {
  
    output_dir_n <- paste(output_dir, p, sep = '_')  

    # Import the node median fluoresencs values and add them together
    table_file_path <- paste(data_file_path, '.density.fcs.cluster.fcs.anno.Rsave_table.csv', sep='')
    table_path <- file.path(output_dir_n,'tables', 'bySample', table_file_path)
    table_path <- normalizePath(table_path)
    SPADE_nodes <- as.data.frame(read.csv(table_path))
    
  if(p == 1) {
    
    combined_SPADE_nodes <- SPADE_nodes
    print(SPADE_nodes[5,9])
    
    for(i in 1:length(SPADE_nodes$ID)){
      
      replicates[[i]] <- list(SPADE_nodes[i,9:14])
      
    }
    
  }
  
  if(p > 1) {
    
    count <- count + 1
    combined_SPADE_nodes <- combined_SPADE_nodes + SPADE_nodes
    print(SPADE_nodes[5,9])
    
    
    for(i in 1:length(SPADE_nodes$ID)){
      
      
      replicates[[i]][[count]] <- SPADE_nodes[i,9:14]
      
    }
  
  }
  }

  print(head(replicates))
  # Calculate the CVs
  
  if(is.list(annotations) == TRUE){
    
    for(i in 1:length(annotations)){
      for(k in 1:length(annotations[[i]])){
      
    }
      
      print(replicates[[annotations[[i]][k]]])
      
    }
    
  }
  
  #average the n combined tables  
  av_SPADE_nodes <- combined_SPADE_nodes / n
  print(n)
  print(combined_SPADE_nodes[5,9])
  print(av_SPADE_nodes[5,9])
  
  # Prune the table values so just the ones we want
  # Use markers from input
  if(panel == 'MSC'){
    cluster_panel <- c('ID', 'mediansComp.APC.A', 'mediansComp.PE.A', 'mediansComp.Pacific.Blue.A', 'mediansComp.APC.Cy7.A', 'mediansComp.FITC.A')
  }
  if(panel == 'TAM'){
    cluster_panel <- c('ID', 'mediansComp.APC.A', 'mediansComp.PE.A', 'mediansComp.Pacific.Blue.A', 'mediansComp.APC.Cy7.A', 'mediansComp.FITC.A', 'mediansComp.AmCyan.A')
    
  }
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
  # So the transform function is dependant on the expression of the marker set
  # This level 0.0095 appears to work alright for separating out the colours
  arcsinh_NODES <- as.data.frame(asinh(prunedByChannel_SPADE_nodes*transform)) 
  
  
  #Rename the columns - this should be done with a double vector of markers and stains
  if(panel == "MSC"){
    colnames(arcsinh_NODES) = c('CD140a', 'Ly6A', 'CD90', 'CD45', 'CD34')
  }
  
  if(panel == 'TAM'){
    colnames(arcsinh_NODES) = c('CD206', 'MHCII', 'F4_80', 'CD11c', 'CD45', 'CD11b')
  }
  
  # Testing print statements
  
  # print(prunedByChannel_SPADE_nodes)
  
  # print(arcsinh_NODES) 
  
  # print(scale_nodes)
  
  # Do the heatmapping
  #Default heat scale- black, yellow, 
  heat_palette <- colorRampPalette(c("black", "yellow", "#FAF7C9"))
  rbPal <- colorRampPalette(c('blue', 'grey', 'red'))
  bl_red <- rbPal(30)
  hot_pal <- heat_palette(80)

  # heat_palette break for ramp
  heat_pairs.breaks <- seq(-1, 7, by=0.1)
  
  # Testing print statements
  print(max(arcsinh_NODES))
  print(min(arcsinh_NODES))
  
  # Red blue heatmap - replace this with the modified breaks to optimize heatmap (min = 0.1, max is 995 percentile or summat)
  blue_breaks <- seq(min(prunedByChannel_SPADE_nodes), max(prunedByChannel_SPADE_nodes), by = (max(prunedByChannel_SPADE_nodes) - min(prunedByChannel_SPADE_nodes)) / 80) 
  blue_breaks <- seq(0 , 8, by = (8 - min(prunedByChannel_SPADE_nodes)) / 80)
  
  # For pheatmap, calculate the optimal fontsize to fit the row
  fontsize_row = 700 / nrow(arcsinh_NODES) 
  
  # Need to fix messed up formatting of dendrogram 
  # shifting to superheat package
  if(is.null(annotations)){
    
    # Genreate unique filename for clustering
    time <- Sys.time()
    time <- format(time, format =  '%H%M')
    pdf_outfile <- paste(output_dir,'asinh', transform, time, 'n', n, '.pdf', sep = '_')

    pdf(pdf_outfile, paper='a4', width = 50, height = 40)


    # pheatmap execution if no annotations
  
    # annotations 
    # annotations <- data.frame(Var1 = factor())
    
    # Create data.frame of annotations
    # Put the annotation function here
    # annotation <- data.frame(Var1 = factor(1:200 %% 2 == 0, labels = c("Exp1", "Exp2")))
  
    # rename the colnames to match the matrix of columns
    # rownames(annotation) <- rownames(arcsinh_NODES)
    
    # Need to create the annotations list to write the clusters onto the heatmap.
    # Booya this works
    
    pheatmap(arcsinh_NODES, 
            color = bl_red, 
            scale = 'none', 
            cellwidth = 8, 
            cellheight =  8,
            fontsize_row = fontsize_row,
            fontsize_col = fontsize_row+1,
            cutree_rows = 10,
            annotation_row = final_annotations)
    
    dev.off()
    
  }
  
  
  # if is annotations then map them back onto the heatmap 
  if (is.list(annotations)) {
      
      # generate vector for categorizing 
      num_cluster <- range(clust_df$cluster)[2]
      
      column_annotations <- seq(1, 80, by=1)

        
        for(k in 1:length(annotations)){
          
          
          # col_annotations[annotations[[k]]] <- k
          positions <- which(column_annotations %in% annotations[[k]])
          col_annotations[positions] <- k
        
          }
        
      
      final_annotations <- as.data.frame(factor(col_annotations, labels = names))
      
      rownames(final_annotations) <- rownames(arcsinh_NODES)
      
      pheatmap(arcsinh_NODES, 
             color = bl_red, 
             scale = 'none', 
             cellwidth = 6, 
             cellheight =  3.5, 
             fontsize_row = fontsize_row,
             row_annotations = col_annotations)
    
    
  }

  
  
  
  
}
