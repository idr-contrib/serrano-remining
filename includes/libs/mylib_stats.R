
# stats.get_FisherFeatureWeights
# stats.get_distanceToControlClassMean
# stats.compute_cosineSimilarityThreshold
# stats.remove_zeroColumns_fromDF
# stats.remove_constantColumns_fromDF
# stats.remove_redundantColumnsByCorrelation
# stats.isConstant_vector
# stats.get_NumberOfClusters_byNullEigenvalues
# stats.get_NumberOfClusters_bySlopeMethod
# stats.get_NearestNeighbors
# stats.testMantel
# stats.summarySE
# stats.calculate_PCA
# stats.calculate_MDS
# stats.calculate_SOM
# stats.calculate_PageRank
# stats.get_normOfAVector
# stats.getdistance_fromPointToSegment
# stats.get_elbow_ofDistribution
# stats.get_numberOfPCAfeaturesToExplain90PercentVariance

############################################################################################
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

stats.remove_constantColumns_fromDF <- function(dfFeatures)
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

stats.get_NumberOfClusters_byNullEigenvalues <- function(laplacian)
{
    # Eigenvalue
    evL <- eigen(L)
    # Determine the number of clusters (null eigenvalues)
    evLRe <- rev(Re(evL$values))
    evLIm <- rev(Im(evL$values))
    less0001Re <- evLRe[evLRe<0.001] 
    less0001Im <- evLIm[evLIm<0.001] 
    numberOfNullEigenvalues <- length(less0001Re) 
    plot(1:1000, rev(Re(evL$values))[1:1000], main=paste0(numberOfNullEigenvalues, " null eigenvalues (",nrow(laplacian)," terms)"), xlab="Eigenvalues")
    abline(v=numberOfNullEigenvalues, col="red", lty=2) 
    return(numberOfNullEigenvalues)
}

############################################################################################

stats.get_NumberOfClusters_bySlopeMethod <- function(eigenvalues)
{
    # If there are no null eigenvalues
    delta <- 2
    tan <- c()
    range <- seq(3,nrow(eigenvalues)-2)
    for (i in range)
    {
      l <- eigenvalues[i,1]
      m1  <- (l - eigenvalues[i-2,1]) / delta
      m2  <- (eigenvalues[i+2,1] - l) / delta 
      tan <- c(tan, abs( (m1-m2) / (1+m1*m2) ))
    }
  
    plot(tan)
    firstMax <- max(tan)
  
    return(firstMax)
}

############################################################################################

stats.get_NearestNeighbors <- function(df, g, measure)
{
    data <- df[c("V1", "V2", measure)]
    # Neighbors of gene g for the given measure
    (neighbors <- subset(data, (data$V1==g | data$V2==g)))
    # Which is the maximum value of measure for these pairs? 
    max <- max(neighbors[measure])
    maximalSet <- neighbors[neighbors[measure]==max,]
    maximalSet[measure] <- 1
    # Return pairs with the maximum values for gene g
    return(maximalSet)
}

############################################################################################

stats.testMantel <- function(m1, m2)
{
    # m1 and m2 are the matrices to be compared
    library(vegan)
    test <- mantel(m1, m2, method="spearman", parallel=3)
    return(test)
}

############################################################################################

stats.summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=F, conf.interval=.95) 
{
    ## Summarizes data.
    ## Gives count, mean, standard deviation, standard error of the mean, and confidence
    ## interval (default 95%).
    ##   data: a data frame.
    ##   measurevar: the name of a column that contains the variable to be summariezed
    ##   groupvars: a vector containing names of columns that contain grouping variables
    ##   na.rm: a boolean that indicates whether to ignore NA's
    ##   conf.interval: the percent range of the confidence interval (default is 95%)
  
    library(doBy)
    
    length2 <- function (x, na.rm=F) 
    {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
  
    # Collapse the data
    formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
    datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
    # Rename columns
    names(datac)[ names(datac) == paste0(measurevar, ".mean") ]    <- measurevar
    names(datac)[ names(datac) == paste0(measurevar, ".sd") ]      <- "sd"
    names(datac)[ names(datac) == paste0(measurevar, ".length2") ] <- "N"
  
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}
  
############################################################################################

stats.calculate_PCA <- function(table)
{
    # PCA (without transformation)
    pca <- prcomp(table, scale.=T, center=T)
    plot(pca, type='l')
    pcaValues <- as.data.frame(pca$x[,1:3])  
    return(pcaValues)
}

############################################################################################

stats.calculate_MDS <- function(distMatrix)
{      
    # Multidimensional scaling
    fit <- cmdscale(distMatrix, eig=T, k=2)
    positions <- as.data.frame(fit$points)
    colnames(positions) <- c("x","y")
    positions$label <- simMatrix[,1]
    return(positions)
}

############################################################################################

stats.calculate_PageRank <- function(simMatrix)
{
    library(igraph)
    g <- graph.adjacency(simMatrix)
    pageRankOutput <- page.rank(g)

    df <- data.frame(pageRankOutput$vector)
    df$node <- names(pageRankOutput$vector)
    return(df)
}

############################################################################################

stats.get_normOfAVector <- function(x1, y1, x2, y2) 
{
  return(sqrt((x2-x1)^2+(y2-y1)^2))
}

############################################################################################

stats.getdistance_fromPointToSegment <- function(test.x, test.y, x1, y1, x2, y2) 
{
  ## Distance from (test.x,test.y) to the segment x1,y1,x2,y2
  distance <- NULL
  intersectingPoint.x <- 0
  intersectingPoint.y <- 0
  
  normOfTheSegment <- stats.get_normOfAVector(x1, y1, x2, y2)
  if( normOfTheSegment < 0.00000001) 
  {
    warning("Short segment")
    return(9999)
  }
  u <- (((test.x - x1) * (x2 - x1)) + ((test.y - y1) * (y2 - y1)))
  u <- u / (normOfTheSegment * normOfTheSegment)
  if((u < 0.00001) || (u > 1)) 
  {
    ## closest point does not fall within the line segment, take the shorter distance to an endpoint
    intersectingPoint.x <- stats.get_normOfAVector(test.x, test.y, x1, y1)
    intersectingPoint.y <- stats.get_normOfAVector(test.x, test.y, x2, y2)
    if (intersectingPoint.x > intersectingPoint.y)
    {
      distance <- intersectingPoint.y
    } else 
    {
      distance <- intersectingPoint.x
    }
  } else 
  {
    ## Intersecting point is on the line, use the formula
    intersectingPoint.x <- x1 + u * (x2 - x1)
    intersectingPoint.y <- y1 + u * (y2 - y1)
    distance <- stats.get_normOfAVector(test.x, test.y, intersectingPoint.x, intersectingPoint.y)
  }
  return(distance)
}

############################################################################################

stats.get_elbow_ofDistribution <- function(xPoints, yPoints)
{
  xyCoordinates <- data.frame(x=xPoints, y=yPoints)
  # Points
  point_upLeft    <- as.vector(xyCoordinates[1,])
  point_downRight <- as.vector(xyCoordinates[nrow(xyCoordinates),])
  # Segments
  diagonalSegment <- point_downRight - point_upLeft
  diagonalSegmentNormalized <- diagonalSegment/sqrt(sum(diagonalSegment^2))
  
  xyCoordinates$distToSegment <- apply(xyCoordinates, 1, function(row) { stats.getdistance_fromPointToSegment(row[1], row[2], point_upLeft[1]$x, point_upLeft[2]$y, point_downRight[1]$x, point_downRight[2]$y)})
  
  maxDistanceToTheSegment <- xyCoordinates[which.max(xyCoordinates$distToSegment),]$x
  return(maxDistanceToTheSegment)
}

############################################################################################

stats.get_numberOfPCAfeaturesToExplain90PercentVariance <- function(pcaObject)
{
  varianceExplained <- summary(pcaObject)$importance[3,]
  numberOfComponentsToExplain90 <- length(varianceExplained[varianceExplained<=0.9])
  return (numberOfComponentsToExplain90)
}

############################################################################################
############################################################################################
############################################################################################


