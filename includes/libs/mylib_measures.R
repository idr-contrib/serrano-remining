
# measures.get_cosineSimilarityMatrix
# measures.get_cosineDistanceMatrix
# measures.get_weightedAngularDistance
# measures.get_weightedAngularSimilarityMatrix
# measures.get_euclideanDistanceMatrix
# measures.get_euclideanSimilarityMatrix
# measures.get_jaccardDistanceMatrix
# measures.get_jaccardSimilarityMatrix
# measures.get_correlationSimilarityMatrix
# measures.get_correlationDistanceMatrix
# measures.get_TFIDFSimilarityMatrix
# measures.get_HammingSimilarityMatrix
# measures.get_cohenKappaSimilarityBetweenTwoVectors
# measures.get_cohenKappaSimilarity_profiles

############################################################################################
############################################################################################

measures.get_cosineSimilarityMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    library(lsa, quietly=T)
    cosineSimilarityMatrix <- cosine(t(data.matrix(matrixWithVariablesToMeasureAsRows)))
    return(cosineSimilarityMatrix)
}

############################################################################################

measures.get_cosineDistanceMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    cosineSimilarityMatrix <- measures.get_cosineSimilarityMatrix(matrixWithVariablesToMeasureAsRows)  
    cosineDistanceMatrix   <- 1 - cosineSimilarityMatrix
    return(cosineDistanceMatrix)
}

############################################################################################

measures.get_weightedAngularDistance <- function (v1, v2, pcaEigenvalues)
{
    num <- sum((v1*v2)/sqrt(pcaEigenvalues))
    den <- sqrt(sum(v1*v1)*sum(v2*v2)) 
    
    return(num/den)
}

############################################################################################

measures.get_weightedAngularSimilarityMatrix <- function (matrixWithVariablesToMeasureAsRows, pcaEigenvalues)
{
    m <- matrixWithVariablesToMeasureAsRows
    allPairs <- t(combn(rownames(m), 2))
    distances <- apply(allPairs, 1, function(x) measures.get_weightedAngularDistance(m[x[1],], m[x[2],], pcaEigenvalues))
    pair_distance <- data.frame(allPairs, distances)
    # Convert to matrix
    source("libs/mylib_utils.R")
    mat <- utils.convert_ListOfPairsToMatrix(pair_distance)
    
    return(mat)
}

############################################################################################

measures.get_euclideanDistanceMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    library(stats)
    return(as.matrix(dist(matrixWithVariablesToMeasureAsRows, method="euclidean")))
}

############################################################################################

measures.get_euclideanSimilarityMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    library(stats)
    return(as.matrix(1 / (1 + dist(matrixWithVariablesToMeasureAsRows, method="euclidean"))))
}

############################################################################################

measures.get_jaccardDistanceMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    require(prabclus, quietly=T)
    jaccardDistanceMatrix <- jaccard(t(matrixWithVariablesToMeasureAsRows))
    return(jaccardDistanceMatrix)
}

############################################################################################

measures.get_jaccardSimilarityMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    jaccardDistanceMatrix <- measures.get_jaccardDistanceMatrix(matrixWithVariablesToMeasureAsRows)
    jaccardSimilarityMatrix <- 1 - jaccardDistanceMatrix
    return(jaccardSimilarityMatrix)
}

############################################################################################

measures.get_correlationSimilarityMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    correlationSimilarityMatrix <- cor(t(matrixWithVariablesToMeasureAsRows))
    return(correlationSimilarityMatrix)
}

############################################################################################

measures.get_correlationDistanceMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    correlationSimilarityMatrix <- measures.get_correlationSimilarityMatrix(matrixWithVariablesToMeasureAsRows)
    correlationDistanceMatrix <- 1 - correlationSimilarityMatrix
    return(correlationDistanceMatrix)
}

############################################################################################

measures.get_TFIDFSimilarityMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
    m <- matrixWithVariablesToMeasureAsRows
    
    idf <- log(nrow(m)/(1+colSums(m)))
    m.idf <- mapply("*", as.data.frame(m), idf)
    rownames(m.idf) <- rownames(m)
    
    (matrix <- t(m.idf))
    nvar <- ncol(matrix)
    simMat <- matrix(nrow=nvar, ncol=nvar)
    for(i in 1:nvar)  for (j in 1:nvar)
    {
        (profile1 <- matrix[,i])
        (profile2 <- matrix[,j])
        intersection <- profile1[profile1==profile2]
        simMat[i,j] <- max(intersection)
    }
    dimnames(simMat) <- list(names(matrix), names(matrix))
    return(simMat)
}

############################################################################################

measures.get_HammingSimilarityMatrix <- function (matrixWithVariablesToMeasureAsRows)
{
	m <- matrixWithVariablesToMeasureAsRows
    nn <- nrow(matrixWithVariablesToMeasureAsRows)
    simMat <- matrix(nrow=nn, ncol=nn)
    for(i in seq_len(nn - 1))
    {
        for(j in seq(i, nn))
        {
            simMat[j,i] <- simMat[i,j] <- sum(m[i,] == m[j,])
        }
    }
    simMat <- simMat/ncol(m)
    return(simMat)
}

############################################################################################

measures.get_cohenKappaSimilarityBetweenTwoVectors <- function(x, y)
{
    # This is a modified version of Cohen's Kappa similarity.
    # Here, 0 in a profile means not tested/not annotated instead of not shown.
    (x0 <- x == 0)
    (x1 <- x == 1)
    (y0 <- y == 0)
    (y1 <- y == 1)
    
    (a <- sum(x1 & y1))
    (b <- sum(x0 & y1))
    (c <- sum(x1 & y0))
    (d <- sum(x0 & y0))
    (N <- sum(a, b, c, d))
    
    # k = (number of 1s in common - number of 1s in common expected by chance) / (profile length - number of 1s in common expected by chance)
    # number of 1s in common by chance = (number of 1s in profile1 * number of 1s in profile2 + number of 0s in profile1 * number of 0s in profile2)/(profile length)^2
    
    numberOf1sInCommon <- a
    numberOf1sInProfile1 <- a+c
    numberOf1sInProfile2 <- a+b
    numberOf0sInProfile1 <- b+d
    numberOf0sInProfile2 <- c+d
    pe <- (numberOf1sInProfile1*numberOf1sInProfile2 + numberOf0sInProfile1*numberOf0sInProfile2)/N^2
    k <- (numberOf1sInCommon-pe) / (N-pe)
    
    return(k)
}

############################################################################################

measures.get_cohenKappaSimilarity_profiles <- function(matrix)
{
    nvar <- ncol(matrix)
    simMat <- matrix(nrow=nvar, ncol=nvar)
    for(i in 1:nvar)  for (j in 1:nvar)
    {
        simMat[i,j] <- measures.get_cohenKappaSimilarityBetweenTwoVectors(matrix[,i], matrix[,j])
    }
    dimnames(simMat) <- list(names(matrix), names(matrix))
    return(simMat)
}

############################################################################################


