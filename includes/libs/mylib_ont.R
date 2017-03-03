
# ont.generate_AnnotationFileFromMatrix

############################################################################################

ont.generate_AnnotationFileFromMatrix <- function(matrix_gene_term)
{
    # Reshape matrix to get a list
    library(reshape2)
    gene_term <- melt(as.matrix(matrix_gene_term))
    # Filter 0 values
    gene_term <- gene_term[gene_term$value!=0,]
    gene_term <- gene_term[c("Var1","Var2")]
    colnames(gene_term) <- c("gene", "term")
    # Split name of the term to get the id, since it initially has the format "name (id:XXXXXXX)"
    termIDS <- data.frame(do.call('rbind', strsplit(as.character(gene_term$term), '(', fixed=T)))
    termIDS$X2 <- substr(termIDS$X2, 1, nchar(as.character(termIDS$X2))-1)
    gene_term$term <- termIDS$X2
    # Generate the list of strings with the annotations
    goaList <- apply(gene_term, 1, function (x) paste("Geneid\t",x[1],"\tXXX00\t\t",x[2],"\tPMID:00000000\tEXP\tGeneid:X00000\tF\tname\tsynonyms\tprotein\ttaxon:9606\t20110627\tIntAct\tGeneid:000-0"))

    return(goaList)
}

############################################################################################



############################################################################################



############################################################################################



############################################################################################
