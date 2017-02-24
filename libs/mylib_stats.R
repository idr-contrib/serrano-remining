############################################################################################

stats.get_FisherFeatureWeights <- function(data, indexID, indexClass, featureRange, LABEL_CONTROL_CLASS)
{
  library(matrixStats, quietly=T)
  
  imageIDS     <- data[indexID][,1]
  featureNames <- colnames(data)[featureRange]
  class        <- data[indexClass]
  # MATRIX (rows=images, cols=features)
  matrix <- data[featureRange]
  rownames(matrix) <- imageIDS
  colnames(matrix) <- featureNames
  # Normalize matrix
  matrix <- scale(matrix, center=F, scale=colSums(matrix))
  submatrixClass <- matrix[class==LABEL_CONTROL_CLASS,]
  num <- (colMeans(matrix)-colMeans(submatrixClass))^2
  den <- colVars(submatrixClass)
  w <- num/den
  
  return(w)
}

############################################################################################

stats.get_distanceToControlClassMean <- function (controlClassFeatures, testClassFeatures)
{
  library(lsa, quietly=T)
  meanControlClass <- apply(controlClassFeatures, 2, mean)
  testImages <- rownames(testClassFeatures)
  distance <- c()
  for (imageIndex in testImages)
  {
    distance <- c(distance, 1 - cosine(as.numeric(testClassFeatures[imageIndex,]), meanControlClass))
  }
  distanceToControlClass <- data.frame(testImages, distance)
  return(distanceToControlClass)
}

############################################################################################

stats.compute_cosineSimilarityThreshold <- function (upperTriangularCosineMatrix)
{
  return(mean(upperTriangularCosineMatrix)-sd(upperTriangularCosineMatrix))
}

############################################################################################

stats.remove_zeroColumns_fromDF <- function(dfFeatures)
{
  features_without0 <- dfFeatures[, which(!apply(dfFeatures, 2, FUN=function(x) {all(x==0)}))]
  return(features_without0)
}

############################################################################################

stats.remove_standardDeviation0_fromDF <- function(dfFeatures)
{
  features_noConstants <- dfFeatures[, which(!apply(dfFeatures, 2, FUN=function(x) {sd(x)==0}))]
  return(features_noConstants)  
}

############################################################################################

stats.remove_redundantColumnsByCorrelation <- function(dfFeatures, cutoffCor)
{
  correlationMatrix <- cor(dfFeatures)
  highlyCorrelatedFeatures <- findCorrelation(correlationMatrix, cutoff=cutoffCor, names=T)
  featuresNamesToBeSelected <- rownames(correlationMatrix)[!(rownames(correlationMatrix) %in% highlyCorrelatedFeatures)]
  return(featuresNamesToBeSelected)
}

############################################################################################

stats.isConstant_vector <- function(v)
{
  return(sd(v)==0)  
}


############################################################################################
