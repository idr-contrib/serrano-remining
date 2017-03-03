
# images.remove_siRNAWithLessThan3ImagesLeft
# images.get_geneSimilarityMatrix_fromAdjacencyMatrixOfImages
# images.get_FeatureVector_fromIndividualFile
# images.add_CharactersToPlateID
# images.add_CharactersToWellID
# images.get_specificPhenotype

############################################################################################
############################################################################################

images.remove_siRNAWithLessThanXImagesLeft <- function (df, X)
{
  freqSiRNA <- table(df$siRNA)
  siRNAselected <- names(freqSiRNA[freqSiRNA>(X-1)])
  df <- subset(df, siRNA %in% siRNAselected)
  return(df)
}

############################################################################################

images.get_geneSimilarityMatrix_fromAdjacencyMatrixOfImages <- function(data, imAdjMatrix)
{
    library(igraph)
    source("libs/mylib_utils.R")
    pairsA <- utils.convert_MatrixToListOfPairs(imAdjMatrix)
    pairsA <- pairsA[pairsA$value!=0, c("X1", "X2")]
    
    (dict <- setNames(data$symbol, data$id))
    
    # Change level to add symbols
    levels(pairsA$X1) <- c(levels(pairsA$X1), levels(dict))
    levels(pairsA$X2) <- c(levels(pairsA$X2), levels(dict))
    
    for (i in 1:nrow(pairsA))
    {
      pairsA$X1[[i]] <- dict[[as.character(pairsA$X1[i])]]
      pairsA$X2[[i]] <- dict[[as.character(pairsA$X2[i])]]
    }
    
    g <- graph.data.frame(pairsA, directed=F)
    #plot(g)
    
    geneSimilarityMatrix <- get.adjacency(g)
    
    return(as.matrix(geneSimilarityMatrix))
}

############################################################################################

images.get_FeatureVector_fromIndividualFile <- function(filepath)
{
    if (file.exists(filepath))
    {
        lineWithFeatures <- 2923
        file <- readLines(filepath)
        ### HEADER
        numberOfClasses <- as.numeric(strsplit(file[1],"\t")[[1]][1]) # 1 class
        numberOfFeatures <- as.numeric(file[2]) # 2919
        numberOfImages <- as.numeric(file[3])   # 1
        featureNames <- read.table(text=file[4:(lineWithFeatures-1)], sep="\n", stringsAsFactors=F)
        if (nrow(featureNames) != numberOfFeatures)
        {
          cat("ERROR: Feature names NOT loaded")
        }
        ### CLASSES
        #classNames <- read.table(text=file[(lineWithFeatures+1):(lineWithFeatures+numberOfClasses)], sep="\n", stringsAsFactors=F)
        ### VECTORS AND IMAGE NAMES
        #featureVectorAndImageNames <- read.table(text=file[(lineWithFeatures+numberOfClasses+1):length(file)], sep="\n", stringsAsFactors=F)
        featureVector <- read.table(text=file[(lineWithFeatures+numberOfClasses+1)], sep="\n", stringsAsFactors=F)
        featureVector <- featureVector$V1
        
        #oddIndex <- 1 #seq(1,nrow(featureVectorAndImageNames),2)
        featureVector <- as.numeric(strsplit(featureVector, " ")[[1]])
        class(featureVector)
        head(featureVector)
        #evenIndex <- 2 #seq(2,nrow(featureVectorAndImageNames),2)
        #imageNames <- as.character(featureVectorAndImageNames[evenIndex,])
        names(featureVector) <- featureNames
        return(featureVector)
    }
}

############################################################################################

images.add_CharactersToPlateID <- function(plateFactor)
{
    return(sprintf("%04d", plateFactor))
}

############################################################################################

images.add_CharactersToWellID <- function(wellFactor)
{
    return(paste0("W", sprintf("%05d", wellFactor)))
}

############################################################################################

images.get_specificPhenotype <- function(ICFilePath, metadata)
{
    # Clean and extract Term-IC associations
    term_IC <- read.csv(ICFilePath, header=F, sep="\t", check.names=F)
    colnames(term_IC) <- c("phenotype", "IC", "tag")
    term_IC <- term_IC[term_IC$IC != "Inf",]
    term_IC$IC <- round(term_IC$IC,3)
    termIDS <- data.frame(do.call('rbind', strsplit(as.character(term_IC$phenotype), '(', fixed=T)))
    termIDS$X2 <- substr(termIDS$X2, 1, nchar(as.character(termIDS$X2))-1)
    term_IC <- data.frame(specificCMPO=termIDS$X2, specificName=termIDS$X1, IC=term_IC$IC, tag=term_IC$tag)
  
    # Split the list of phenotypes
    metadata$splitCMPO <- lapply(as.character(metadata$`CMPO ID(s)`), function(x) strsplit(x, "|", fixed=T))
    
    # Get the list siRNA-name of the individual phenotype
    siRNA_term <- subset(metadata, select=c(siRNA, `CMPO term(s)`))
    splitPhenotypes <- strsplit(as.character(metadata$`CMPO term(s)`), "|", fixed=T)
    edgelist <- unique(data.frame(siRNA=rep(siRNA_term$siRNA, sapply(splitPhenotypes,length)), ph=unlist(splitPhenotypes)))
  
    # Most specific term for each siRNA
    specificIDs <- c()
    specificNames <- c()
    for (phList in metadata$splitCMPO)
    {
        df <- data.frame(ph=phList)
        colnames(df) <- "specificCMPO"
        df <- unique(merge(term_IC, df))
        if (nrow(df)==0)
        {
            term <- "nophenotype"
            name <- "nophenotype"
        }
        else
        {
            # Term with max IC
            term <- as.character(df[df$IC==max(df$IC),]$specificCMPO)
            name <- as.character(df[df$specificCMPO==term,]$specificName)
        }
        specificIDs   <- c(specificIDs, term)
        specificNames <- c(specificNames, name)
    }
    
    metadata$specificCMPO <- specificIDs
    metadata$specificName <- specificNames
    
    # Add the IC and tag for each phenotype
    metadata <- merge(metadata, term_IC, all.x=T)
    
    return(metadata)
}

############################################################################################


