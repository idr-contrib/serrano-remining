
# utils.convert_MatrixToListOfPairs
# utils.convert_ListOfPairsToMatrix
# utils.remove_NA_Inf_0_fromVector
# utils.remove_elementsLowerThanThreshold_fromVector
# utils.get_highestValues_fromAVector
# utils.scale_featureUsingWeights_fromDF
# utils.sort_vector
# utils.split_stringList_to_df
# utils.extract_IDsFromNameList
# utils.normalize
# utils.mappingENSGtoEntrez_Biomart
# utils.vector2string
# utils.differenceOfDataFrames

############################################################################################
############################################################################################

utils.convert_MatrixToListOfPairs <- function(matrix)
{
    library(reshape, quietly=T)
    
    matrix[!lower.tri(matrix)] <- NA
    pairList <- melt(matrix, na.rm=TRUE)
    pairList <- pairList[!is.na(pairList$value),]
    return(pairList)
}

############################################################################################

utils.convert_ListOfPairsToMatrix <- function(df_pair_value)
{
    library(igraph)
  
    g <- graph.data.frame(df_pair_value)
    E(g)$weight <- df_pair_value[,3]
    mat <- get.adjacency(g, sparse=FALSE, attr="weight")
    return(mat)
}

############################################################################################

utils.remove_NA_Inf_0_fromVector <- function(vector)
{
    vector <- vector[!is.na(vector)]
    vector <- vector[vector!="Inf"]
    vector <- vector[vector!=0]
    return(vector)
}

############################################################################################

utils.remove_elementsLowerThanThreshold_fromVector <- function(vector, threshold)
{
    vector <- vector[vector>=threshold]
    return(vector)
}

############################################################################################

utils.get_highestValues_fromAVector <- function(vector, percentageTopValues)
{
    vector  <- vector[order(vector, decreasing=T)]
    filteredVector <- head(vector, percentageTopValues*length(vector))
    return(filteredVector)
}

############################################################################################

utils.scale_featureUsingWeights_fromDF <- function(dfFeatures, weights)
{
    weightedFeatures <- t(t(dfFeatures) * weights)
    return(weightedFeatures)
}

############################################################################################

utils.sort_vector <- function(v)
{
    return(v[order(v)])
}

############################################################################################

utils.split_stringList_to_df <- function(str, delim)
{
    return(data.frame(do.call('rbind', strsplit(str,delim,fixed=T))))
}

############################################################################################

utils.extract_IDsFromNameList <- function(listWithNameAndID)
{
    termIDS <- data.frame(do.call('rbind', strsplit(as.character(listWithNameAndID), '(', fixed=T)))
    code <- substr(termIDS$X2, 1, nchar(as.character(termIDS$X2))-1)
    return(code)
}

############################################################################################

utils.normalize <- function(vector)
{
    return( (vector - min(vector))/(max(vector)-min(vector)) )
}

############################################################################################

utils.mappingENSGtoEntrez_Biomart <- function(ENSG_ids)
{
    library("biomaRt")
    mart <- useMart("ensembl")
    mart <- useDataset("hsapiens_gene_ensembl", mart)
    ensembl2entrez  <- getBM(attributes=c("ensembl_gene_id", "entrezgene"), filters="ensembl_gene_id", values=ENSG_ids, mart=mart)
    ensembl2entrez <- ensembl2entrez[!is.na(ensembl2entrez$entrezgene),]
    return(unique(ensembl2entrez))
}

############################################################################################

utils.vector2string <- function(vector)
{ 
    return(paste(vector, sep="", collapse=""))
}

############################################################################################

utils.differenceOfDataFrames <- function(allData, dfToRemove)
{
    require(dplyr) 
    cleanData <- anti_join(allData, dfToRemove)
    return(cleanData)
}

############################################################################################

############################################################################################
