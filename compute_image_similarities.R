library(optparse)
library(rhdf5)
library(proxy)


option_list = list(
  make_option(c("-i", "--input"), action = "store", default = NA, type='character',
              help="Path to the HDF5 file with the WND-CHRM features"),
  make_option(c("-o", "--output"), action = "store", default = NA, type='character',
              help="Basename to use for the output files (PCs and similarity matrix)"),
  make_option(c("-n", "--nPCs"), action = "store", default = NA, type = "integer",
              help = "Number of principal components to keep and use")
)

opt_parser <- OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.na(opt$input)){
  print_help(opt_parser)
  stop("Argument --input required.", call.=FALSE)
}


# Helper function to find elbow in a 2D curve represented by a list of ordered values
# This function finds the point with maximum distance from the line between
# the first and last points.
# Adapted from StackOverflow:
# http://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve

find_elbow <- function(values) {
  
  n <- length(values)
  
  # Crude check to see if values are ordered
  if (values[1]<=values[n]) {
    stop("Values must be in decreasing order") 
  }
  
  coords <- cbind(seq.int(n), as.matrix(values))
  # Vector between first and last point
  line <- coords[n,] - coords[1,]
  # Normalize the vector
  line <- line/sqrt(sum(line^2));
  # Vector between all points and first point
  V1 <- sweep(coords, 2, coords[1,]) # default function in sweep is - (minus)
  # To calculate the distance to the line, we split V1 into two
  # components, V2 that is parallel to the line and V3 that is perpendicular.
  # The distance is the norm of the part that is perpendicular to the line.
  # We find V2 by projecting V1 onto the line by taking the scalar product
  # of V1 with the unit vector on the line (this gives us the length of the
  # projection of V1 onto the line).
  # We get V2 by multiplying this scalar product by the unit vector.
  # The perpendicular vector is V1 - V2
  V2 <- (V1 %*% line) %*% line
  V3 <- V1 - V2
  # Distance to line is the norm of V3
  dist <- sqrt(rowSums(V3^2))
  idx <- which.max(dist)
  
  return(coords[idx,])
}


# Input file with WND_CHRM features
filename <- opt$input
fields <- h5ls(filename)

# Read HDF5 file
group_name <- paste(fields$group, fields$name, sep ="/")
h5data <- h5read(filename, group_name[1], compoundAsDataFrame=FALSE)
H5close()

# Get metadata
measures <- h5data$Measurements
imageID <- measures$ImageID
wellID <- measures$WellID

# Features are stored in a list of matrices
listOfFeatureMatrices <- measures[11:length(measures)]

# Make matrix of images x features
featureMatrix <- t(do.call(rbind, listOfFeatureMatrices))
row.names(featureMatrix) <- imageID

# There are multiple rows with the same imageID because of tiling
# We average all rows from the same image
featureMatrix <- aggregate(featureMatrix, by = list(row.names(featureMatrix)), mean)
rownames(featureMatrix) <- featureMatrix[,1]
featureMatrix <- featureMatrix[,-1]

# Remove constant features
featureMatrix <- featureMatrix[, which(!apply(featureMatrix, 2, FUN=function(x) {sd(x)==0}))]
rm(fields, group_name, h5data, measures, imageID, wellID, listOfFeatureMatrices)

# PCA
pca <- prcomp(featureMatrix, scale.= TRUE, center = TRUE)

# Select PCs
k <- NULL
if (!is.na(opt$nPCs)) {
  k <- opt$nPCs
} else {
  k <- find_elbow(pca$sdev)[1]
}
saveRDS(pca$x[,1:k], paste0(opt$output,"-PCs.rds"))

# Cosine similarity
S <- as.matrix(simil(pca$x[,1:k], method = "cosine", by_rows = TRUE), diag = 1)
saveRDS(S, paste0(opt$output,"-similarities.rds"))
