library(SPADE)


# Import marker panels via a list values

# Example, it appears that these should be from the FCS, specify which markers are to be used for clustering
# markers <- c("Cd(110,111,112,114)","Cell_length","Dy(163.929)-Dual","Er(165.930)-Dual","Er(166.932)-Dual","Er(167.932)-Dual","Er(169.935)-Dual","Eu(150.919)-Dual","Eu(152.921)-Dual","Gd(155.922)-Dual","Gd(157.924)-Dual","Gd(159.927)-Dual","Ho(164.930)-Dual","In(114.903)-Dual","Ir(190.960)-Dual","La(138.906)-Dual","Lu(174.940)-Dual","Nd(141.907)-Dual","Nd(143.910)-Dual","Nd(144.912)-Dual","Nd(145.913)-Dual","Nd(147.916)-Dual","Nd(149.920)-Dual","Pr(140.907)-Dual","Sm(146.914)-Dual","Sm(151.919)-Dual","Sm(153.922)-Dual","Tb(158.925)-Dual","Tm(168.934)-Dual","Yb(170.936)-Dual","Yb(171.936)-Dual","Yb(173.938)-Dual","Yb(175.942)-Dual")

marker <- c()

# Specify the 'PANELS'; i.e the sample files (FCS) upon which the marker info will be extracted and run
# Example
# PANELS <- list(list(panel_files=c("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs"), median_cols=NULL,reference_files=c("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs"),fold_cols=c()))

PANELS <- list(list(panel_files=c(''), median_cols=NULL, reference_files=c(''),fold_cols=c()))

# Run SPADE data table and tree generation via SPADE.driver
# What is the layout function here for?

# SPADE.driver("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs", out_dir="output", cluster_cols=markers, panels=PANELS, transforms=flowCore::arcsinhTransform(a=0, b=0.2), layout=SPADE.layout.arch, downsampling_target_percent=0.1, downsampling_target_number=NULL, downsampling_target_pctile=NULL, downsampling_exclude_pctile=0.01, k=200, clustering_samples=50000)

# Specify the table and fucntion to map the force directed treem this can be either a matrix 
# It appears that layout.kamada.kawai force directed tree work quite well

# layout <- read.table(file.path("output","layout.table"))

# Input the igraph object for mapping into  
# mst <- read.graph(file.path("output","mst.gml"),format="gml")

# Plot the trees for each marker expression density
# Layout at this point seems to specify the way that the force directed graph works

# SPADE.plot.trees(mst,"output",file_pattern="*fcs*Rsave", layout=layout.kamada.kawai ,out_dir=file.path("output","pdf1"),size_scale_factor=1.2)


