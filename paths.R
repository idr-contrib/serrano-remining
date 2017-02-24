
# Path to images
pathTo_source_secretionImages <- file.path("/g", "data", "bio-it_centres_course", "data", "VSVG")

pathTo_root <- "screens"

# First level paths
pathTo_screen <- file.path(pathTo_root, screen)
pathTo_data   <- file.path(pathTo_root, screen, "data")
pathTo_plots  <- file.path(pathTo_root, screen, "plots")

    # Second level paths to data
    pathTo_cosine_images <- file.path(pathTo_data, "cosineImages")
    pathTo_cosine_genes  <- file.path(pathTo_data, "cosineGenes")
    pathTo_metadata      <- file.path(pathTo_data, "metadata")
    pathTo_Symlinks      <- file.path(pathTo_data, "symlinks")
    pathTo_averages      <- file.path(pathTo_data, "averages")
    pathTo_wndchrm       <- file.path(pathTo_data, "wndchrm")
    pathTo_clusteringCC  <- file.path(pathTo_data, "clusteringCC")

        # Third level
        pathTo_featureWeights  <- file.path(pathTo_wndchrm, "featureWeights")

    # Second level paths to plots
    pathTo_controlClass <- file.path(pathTo_plots, "controlClass")
    pathTo_hclust       <- file.path(pathTo_plots, "hclust")
    pathTo_MDS          <- file.path(pathTo_plots, "MDS")


