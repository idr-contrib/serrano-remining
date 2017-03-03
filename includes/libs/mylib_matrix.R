
# matrix.get_upperTriangularPart_ofMatrixValues
# matrix.get_adjacencyMatrixGivenTheThreshold
# matrix.get_sharedNearestNeighborsGraph_fromBinaryMatrix
# matrix.remove_diagonal
# matrix.calculate_simpleLaplacian
# matrix.convert_binaryMatrixToSparseBinaryMatrix

############################################################################################
############################################################################################

matrix.get_upperTriangularPart_ofMatrixValues <- function(matrix)
{
    return(matrix[upper.tri(matrix, diag=F)])  
}

############################################################################################

matrix.get_adjacencyMatrixGivenTheThreshold <- function(matrix, threshold)
{
    return(1.0 * (matrix > threshold))
}

############################################################################################

matrix.get_sharedNearestNeighborsGraph_fromBinaryMatrix <- function(BinaryAdjacencyMatrix)
{
    return(BinaryAdjacencyMatrix %*% BinaryAdjacencyMatrix)
}

############################################################################################

matrix.remove_diagonal <- function(matrix)
{
    return(matrix - diag(nrow(matrix)))
}

############################################################################################

matrix.calculate_simpleLaplacian <- function(AdjacencyMatrix)
{
    # Simple Laplacian: L = I − D^(−1) * AdjacencyMatrix
    # Degree inverse
    Dinv <- diag(apply(AdjacencyMatrix, 1, function(x) 1/sum(x)))
  
    # Normalized Adjacency matrix
    normalizedA <- Dinv %*% AdjacencyMatrix
  
    ### Laplacian
    L <- diag(nrow(AdjacencyMatrix)) - normalizedA # simple Laplacian
    L <- round(L,3)
    colnames(L) <- colnames(AdjacencyMatrix)
    rownames(L) <- rownames(AdjacencyMatrix)
  
    return(L)
}

############################################################################################

matrix.convert_binaryMatrixToSparseBinaryMatrix <- function(m, percRemovingOnes)
{
    numberOfOnes <- sum(m)
    onesToRemove <- round(percRemovingOnes * numberOfOnes)
  
    # get positions with 1s
    indexes <- which(m==1, arr.ind=T)
  
    # change 1 by 0 in onesToRemove random positions
    randomIndexes <- indexes[sample(nrow(indexes)),]
    positionsToRemove <- head(randomIndexes, onesToRemove)
  
    m[positionsToRemove] <- 0

    # Clean matrix in case it has empty rows/columns
    m <- m[!colSums(m)==0]
    m <- m[!rowSums(m)==0,]
    return(m)
}

############################################################################################

############################################################################################

############################################################################################
