hill.div <- function(abund,qvalue,tree){ 
  
#Quality-check and warnings
if(missing(abund)) stop("Abundance data is missing")
if(missing(qvalue)) stop("q value is missing")
if(qvalue==1){qvalue=0.99999} 
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")

#Diversity
if(missing(tree)){
    #If input data is a vector
    if(is.null(dim(abund)) == TRUE){
    if(sum(abund) != "1") stop("The abundance data does not sum up to 1")
    pi <- abund[abund!=0] 
    div <- sum(pi^qvalue)^(1/(1-qvalue))
    return(div) 
    }

    #If input data is an OTU table
    if(is.null(dim(abund)) == FALSE){
    samples <- colnames(abund)

    sample.vector <- c()
    for (s in samples){
      vector <- abund[,s]
      pi <- vector[vector!=0] 
      div <- sum(pi^qvalue)^(1/(1-qvalue))
      names(div) <- s
      sample.vector <- c(sample.vector,div)
    }
    return(sample.vector) 
    }
#Phylodiversity
}else{
    if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")  
      
    #If input data is a vector
    if(is.null(dim(abund)) == TRUE){
    if(identical(sort(names(abund)),sort(tree$tip.label)) == FALSE) stop("OTU names in the vector and tree do not match")  
    if(sum(abund) != "1") stop("The abundance data does not sum up to 1")
    Li <- tree$edge.length #Get branch lengths
    ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node)) #Sum relative abundances per lineage
    ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector]))) #Sum relative abundances per lineage
    T <- sum(Li * ai) #Get total tree depth
    Li <- Li[ai != 0] #Remove zeros
    ai <- ai[ai != 0] #Remove zeros
    phylodiv <- sum(Li/T * ai^qvalue)^(1/(1-qvalue)) #Compute phylodiversity
    return(phylodiv)
    }

    #If input data is an OTU table
    if(is.null(dim(abund)) == FALSE){    
    if(identical(sort(rownames(abund)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")  
    samples <- colnames(abund)
    sample.vector <- c()
    for (s in samples){
        vector <- abund[,s]
        names(vector) <- rownames(abund)
        Li <- tree$edge.length #Get branch lengths
        ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node)) #Sum relative abundances per lineage
        ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector]))) #Sum relative abundances per lineage
        T <- sum(Li * ai) #Get total tree depth
        Li <- Li[ai != 0] #Remove zeros
        ai <- ai[ai != 0] #Remove zeros
        phylodiv <- sum(Li/T * ai^qvalue)^(1/(1-qvalue)) #Compute phylodiversity
        names(phylodiv) <- s
        sample.vector <- c(sample.vector,phylodiv)
        }      
    return(sample.vector) 
    }
                       
}
  
}
