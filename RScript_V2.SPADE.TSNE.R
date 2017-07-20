# JO working script with updates from Shahram
# 01/03/2017
# V2 - integrated functional script

# To Do list:
# Make the heatmap all fancy 
# Make the colour palette really nice
# Re-base 
# Driver funtion

# Need to trace modify the write.FCS fucntion to alter the $PnR - done
# line 35 -> change 1025

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
  
  clust_df <- as.data.frame(tSNE_clust)
  
  
  # define color schema
  # must be same colour length as the number clusters used for SPADE
  # Not this is in a 12 colour rainbow, must figure out how to get mroe colours
  clust_df$Col <- palette(rainbow(200))[as.numeric(cut(clust_df$cluster, breaks = 200))]
  
  time <- Sys.time()
  time <- format(time, format =  '%H%M')
  pdf_outfile <- paste(output_dir,'all_cluster', time, '.pdf', sep = '_')
  print(pdf_outfile)
  pdf(pdf_outfile, width = 12, height = 12)
  
  #plot as before
  plot.default(clust_df$tSNE_X_P_65_E_200_I_1500_T_0.3, 
               clust_df$tSNE_Y_P_65_E_200_I_1500_T_0.3, 
               pch=19, 
               cex=0.3,
               col = clust_df$Col,
               xlab = 'tSNE 1',
               ylab = 'tSNE 2'
               )
  dev.off()
  
  
  # return(clust_df)
  
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

# Correctly annotate the tSNE dotplot with the lcuster annotations
# In progress - non functional
annotate.cluster.plot <- function(cluster_fcs_file, annotations, names){
  
  # Take the heatmap clusters and split them adding annotation
  
  
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
  
  xyplot(Comp.PE.A~Comp.APC.A,ff)
  
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
  library(RColorBrewer)
  
}



QC_fcs_file <- function(data_file_path) {
  
  
  # QC with user interface to check
  fsc_file <- flowCore::read.FCS(data_file_path, truncate_max_range = FALSE)
  fsc_file

}

# This function should make dotplots of the tSNE graphs
# tSNE graphs should be heat cloloured by channel
# With scale transform specified use asinh with t factor = 0.0095

color.bar <- function(lut, pdf_file , min, max=-min, nticks=10, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  scale = round(scale, 3)  
  
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

# This is a function using min/max normalization
# Except going to use normalization to 99%-ile to prune for outliers
# This function can be applited to a df using lapply
normalize <- function(x) {
  return ((x - min(x)) / (quantile(x, 0.995) - min(x)))
}

asinh_scale <- function(x, transform){
  return(asinh(x*transform))
}
# Need to make the transform a vector and apply each transformation to each channel individually
plot_tSNE <- function(csv_file_path, transform) {
  
  # read the accense data in put and plot the flowJO tSNE points
  tSNE_data <- read.csv(csv_file_path, sep = ',')
  
  # Do i need this line of code?
  transform_tSNE <- as.matrix(tSNE_data)
  
  # Create red/blue color ramp
  # spectralPal <- colorRampPalette(brewer.pal(11,"Spectral"))(100)
  # rbPal <- colorRampPalette(c('blue', 'grey', 'red'))
  
  
  count <- 1
  
  # Create a new folder for to put the test in named by the transfrom factor
  fp <- getwd()
  folder <- paste(fp, transform[1], sep = '/')
  
  dir.create(folder)
  
  # Create heat colour scale and add it as a vector column
  x <- ncol(tSNE_data-1)
  y <- ncol(tSNE_data)
  
  for (k in 3:(ncol(tSNE_data)-2)) {
  
    tSNE <- as.data.frame(tSNE_data[k])
    
    # Code to use a vector of data to transform 
    # if (transform != 0){
    
    # lApply arcsinh transformation function
    transform_tSNE <- as.data.frame(lapply(tSNE, asinh_scale, transform[count]))  
    # Normalize the data set with lapply
    transform_tSNE_Norm <- as.data.frame(lapply(transform_tSNE, normalize))
    
    # create a vecor matrix with  
    tran_tSNE_Norm_mt <- as.matrix(transform_tSNE_Norm)
    
    # Assign data to colour bis
    tSNE_cols <- as.numeric(cut(tran_tSNE_Norm_mt, breaks = c(seq(0 , 1, length.out = 100))))
    # Old assignment method
    # tSNE_cols <- as.numeric(cut(transform_tSNE, breaks = c(round(min(transform_tSNE), 2), seq( -0.1, round(quantile(transform_tSNE, 0.9995, na.rm = TRUE), 2), len = 99))))
    # tSNE_cols_2 <- as.numeric(cut(transform_tSNE, breaks = c(round(min(transform_tSNE), 2), seq( -0.1, round(max(transform_tSNE), 2), len = 99))))
    
    # generate matlablike colours
    plot_Cols <- colorRamps::matlab.like2(100)[tSNE_cols]
    
    # Check if all the points are present
    print(length(tSNE_cols))
    print(count)
    print(transform[count])
    
    # Create the output_file 
    pdf_outfile <- paste(colnames(tSNE_data[k]), 'asinh', transform[count], '.pdf', sep = '_')
    path <- paste(folder, pdf_outfile, sep = '/')
    pdf(path, width = 10, height = 10)
    
    count <- count + 1
    # split the pdf into two bits for legend
    # layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
    
    
    # Plot the tSNE
    # Need to change the tSNE parameter name to fi the run parameters
    plot.default(tSNE_data$tSNE_X_P_65_E_200_I_1500_T_0.3, 
                tSNE_data$tSNE_Y_P_65_E_200_I_1500_T_0.3, 
                pch=19, 
                cex=0.5,
                col = plot_Cols,
                xlab = 'tSNE 1',
                ylab = 'tSNE 2'
                )
    
    dev.off()
    
    pdf_bar <- paste(colnames(tSNE_data[k]), 'Legend', '.pdf', sep = '_')
    pdf_outpath <- paste(folder, pdf_bar, sep = '/')
    
    # Create the legend bar
    color.bar(colorRamps::matlab.like2(99), pdf_outpath, min = 0, max = 1)
  
}

}

# Identify number if cell clusters, not 100% sure how many is optimal, 
# Detect number of clusters with NBclust?
# it appears that overclustring is important
# Clustering should be proportional to the cell number as you dont want loads of nodes with no cells in it 
# Does SPADE have to be run with the same random seed everytime for conistency?

RUN_SPADE_tSNE <- function(data_file_path, clusters=100, n=3, transforms='none', markers = c('tSNE_X_P_65_E_200_I_1500_T_0.3', 'tSNE_Y_P_65_E_200_I_1500_T_0.3'))  {
  

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


# Note that names always has to begin with 'Un' if not all subsets are specified
# How does the filepath work?

av_table_heatmap_SPADE <- function(data_file_path, n=1, transform, annotations = NULL, fcs_file_cluster = NULL, names = NULL, panel, cutree = 1 ) {

  # Function to average the three tables
  # Then write out the average table and heatmap it

  out <- strsplit(data_file_path,' ')
  out_path <- out[[1]][1]
  output_dir <- file.path(out_path)  
  replicates <- list()
  count <- 1
  
  for(p in 1:n) {
    print(p)
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
    cluster_panel <- c('ID', 'mediansComp.APC.A', 'mediansComp.PE.A', 'mediansComp.Pacific.Blue.A', 'mediansComp.FITC.A')
  }
  if(panel == 'TAM'){
    cluster_panel <- c('ID', 'mediansComp.APC.A', 'mediansComp.PE.A', 'mediansComp.Pacific.Blue.A', 'mediansComp.APC.Cy7.A', 'mediansComp.AmCyan.A')
    
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
  prunedByChannel_SPADE_nodes <- as.data.frame(prunedByChannel_SPADE_nodes[,2:x])
  
  arcsinh_NODES <- prunedByChannel_SPADE_nodes
  
  count <- 1
  # Transform needs to be a vector
  # Normalize each row here for the heatmap table
  for(k in 1:ncol(prunedByChannel_SPADE_nodes)){
    print(head(prunedByChannel_SPADE_nodes[k]))
    tSNE <- as.data.frame(prunedByChannel_SPADE_nodes[k])
    
    # lApply arcsinh transformation function
    transform_tSNE <- as.data.frame(lapply(tSNE, asinh_scale, transform[count]))  
    
    # Normalize the data set with lapply
    # transform_tSNE_Norm <- as.data.frame(lapply(transform_tSNE, normalize))
    
    arcsinh_NODES[k] <- transform_tSNE    
  
  }
  
  print(head(arcsinh_NODES))
  
  # R Colourbrewer set3
  #Rename the columns - this should be done with a double vector of markers and stains
  if(panel == "MSC"){
    colnames(arcsinh_NODES) = c('CD140a', 'CD31', 'CD90', 'CD34')
    # bl_red_Brown_pink
    my_cols <- c('#377eb8', '#e41a1c', '#d95f02', '#7570b3')
  }
  
  # R colourbrewer set1
  if(panel == 'TAM'){
    colnames(arcsinh_NODES) = c('CD206', 'MHCII', 'F4_80', 'CD11c', 'CD11b')
    # bl_orange_green_turq_cyan
    my_cols <- c('#377eb8')
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
  # heat_pairs.breaks <- seq(-1, 7, by=0.1)
  
  # Testing print statements
  print(max(arcsinh_NODES))
  print(min(arcsinh_NODES))
  
  # Red blue heatmap - replace this with the modified breaks to optimize heatmap (min = 0.1, max is 995 percentile or summat)
  # Create a list of vectors in order to bin the values into for the breaks
  blue_breaks <- seq(min(prunedByChannel_SPADE_nodes), max(prunedByChannel_SPADE_nodes), by = (max(prunedByChannel_SPADE_nodes) - min(prunedByChannel_SPADE_nodes)) / 80) 
  
  # blue_breaks is throwing some error here
  # blue_breaks <- seq(0 , 8, by = (8 - min(prunedByChannel_SPADE_nodes)) / 60)
  
  # For pheatmap, calculate the optimal fontsize to fit the row
  fontsize_row = 700 / nrow(arcsinh_NODES) 
  
  # Need to fix messed up formatting of dendrogram 
  # shifting to superheat package
  if(is.null(annotations)){
    
    # Genreate unique filename for clustering
    time <- Sys.time()
    time <- format(time, format =  '%H%M')
    pdf_outfile <- paste(output_dir,'asinh', 'vector_transform', time, 'n', n, '.pdf', sep = '_')

    pdf(pdf_outfile, width = 12, height = 12)


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
            cellwidth = 15, 
            cellheight =  15,
            fontsize_row = fontsize_row,
            fontsize_col = fontsize_row+0.5,
            cutree_rows = cutree)
            # annotation_row = final_annotations
            
    
    dev.off()
    
  }
  
  
  # if is annotations then map them back onto the heatmap 
  if (is.list(annotations)) {
      
    print('graphing with annotations')
    
    # Find the cluster file in the file
    cluster_file <- fcs_file_cluster
    
    print(cluster_file)
    
    # read and display the cluater_fcs file
    cluster_fcs <- flowCore::read.FCS(cluster_file, truncate_max_range = FALSE)
    
    # extrat the data from within
    tSNE_clust <- exprs(cluster_fcs)
    
    clust_df <- as.data.frame(tSNE_clust)
    
    # generate vector for categorizing 
    num_cluster <- range(clust_df$cluster)[2]
  
    # Do the heatmap annotations  
    col_annotations <- rep(0, num_cluster)

        
      for(k in 1:length(annotations)){
          subset <- annotations[[k]]
          for(m in 1:length(subset)){
            col_annotations[annotations[[k]][m]] <- k 
            
        
          }
      }
      
      final_annotations <- as.data.frame(factor(col_annotations, labels = names))
      
      rownames(final_annotations) <- rownames(arcsinh_NODES)
      
      colnames(final_annotations) <- c('Cell subset')
      
      ann_colours <- list('red', 'blue', 'green', 'yellow')
      
      pdf('annotations.pdf', width = 12, height = 12)
      
      if(panel == 'TAM'){
        palette(brewer.pal(n = 8, name = 'Set1'))
      }
      
      if(panel == 'MSC'){
        palette(brewer.pal(n = 8, name = 'Dark2'))
      }
      
        
      pheatmap(arcsinh_NODES, 
             color = bl_red, 
             scale = 'none', 
             cellwidth = 15, 
             cellheight =  15, 
             cutree_rows = cutree,
             fontsize_row = fontsize_row,
             fontsize_col = fontsize_row+0.25,
             annotation_row = final_annotations,
             annotation_colors = ann_colours)
    
      dev.off()
      
      # Create new vector for cluster ID assignment
      num_cells <- length(clust_df$cluster)
      clust_group <- rep(0, num_cells)
      
      # Re-label the vector with grouping assignments
      for(k in 1:length(annotations)){
        # Pick out each list
        subset <- annotations[[k]]
        
        for(m in 1:length(subset)){
          cells <- which(clust_df$cluster == subset[m], arr.ind = TRUE)
          
          for(x in 1:length(cells)){
            
            clust_group[cells[x]] <- k
          }
          
          
        }
      }
      
      final_clust_annotations <- factor(clust_group, labels = names)
      
      # Turn these assignments into colours 
      cc <- palette()
      clust_df$Col <- cc[(final_clust_annotations)]
      
      print(head(clust_df$Col))
      
      # Plot the cluster on the tSNE
      pdf('tSNE_cluster_assignment.pdf', width = 12, height = 12)
      
      
      #plot as before
      plot.default(clust_df$tSNE_X_P_65_E_200_I_1500_T_0.3, 
                   clust_df$tSNE_Y_P_65_E_200_I_1500_T_0.3, 
                   pch=19, 
                   cex=0.5,
                   col = clust_df$Col,
                   xlab = 'tSNE 1',
                   ylab = 'tSNE 2'
                   )
      dev.off()
    
  }
  
  
  
  
}
