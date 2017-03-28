args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One argument required (h5-file)", call.=FALSE)
}

#### 1) Install and load *rhdf5*
if (!("rhdf5" %in% rownames(installed.packages())))
{
  source("http://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
}

library(rhdf5)

#### 2) Explore the structure
filename <- args[1]
fields <- h5ls(filename)
str(fields)

#### 3) Read *example.h5*
group_name <- paste0(fields$group, fields$name)
data <- h5read(filename, group_name[1], compoundAsDataFrame=FALSE)
H5close()

#### 4) Get metadata
measures <- data$Measurements
imageID <- measures$ImageID
wellID <- measures$WellID

#### 5) Get feature values (imageID: `r imageID`; well: `r wellID`)
# Features are stored in a list of matrices
featureListOfMatrices <- measures[11:length(measures)]

## 1) Feature values
# Feature Vector has 2919 values
featureVector <- as.vector(do.call(rbind, featureListOfMatrices))

## 2) Feature names: Build new ID for each feature
# Length of each feature type
feat_size <- lapply(featureListOfMatrices, length)
feat_size <- data.frame(name=names(feat_size), size=unlist(feat_size))
rownames(feat_size) <- seq(1:nrow(feat_size))
# Repeat "length" times the name
a <- rep(feat_size$name, feat_size$size)
# Build sequences of number to create the ids
b <- c()
for(nElem in feat_size$size)
{
  b <- c(b, seq(1:nElem))
}

(feature_value <- data.frame(feature=paste(a,b,sep="_"), value=featureVector))
