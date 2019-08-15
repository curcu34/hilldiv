#' Hill numbers computation
#' @title Hill numbers computation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV diversity
#' @description Compute neutral or phylogenetic Hill numbers from a single sample (vector) or count table (matrix). Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs, i.e. the number of equally abundant OTUs that would be needed to give the same value of diversity.
#' @param abund A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @param qvalue A positive integer or decimal number (>=0), usually between 0 and 3.
#' @param tree An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).
#' @param type Either "abundance" for abundance-based Hill number computation or "incidence" for incidence-based Hill number computation. In the latter case, the 'abund' object needs to be an OTU table, as incidence-based Hill number computation requires multiple samples. Default: type="abundance".
#' @seealso \code{\link{index.div}}, \code{\link{div.part}}
#' @examples
#' hill.div(otu.vector,0)
#' hill.div(otu.table,1)
#' hill.div(otu.table,1,tree)
#' hill.div(otu.table,1,tree,type="incidence")
#' hill.div(otu.table,qvalue=2,tree=tree)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363–375.\cr\cr
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its con‐ sequences. Ecology, 54, 427–432.
#' @export

hill.div <- function(abund,qvalue,tree,type){ 
  
#Quality-check and warnings
if(missing(abund)) stop("Abundance data is missing")
if(missing(qvalue)) stop("q value is missing")
if(qvalue==1){qvalue=0.99999} 
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(type)){type="abundance"}
if((type=="incidence") & (length(dim(abund)) != 2)) stop("Incidence-based Hill number computation requires multiple samples")
  
#Diversity
if(missing(tree)){
    #If input data is a vector
    if(is.null(dim(abund)) == TRUE){
    if(type == "incidence") stop("Incidence-based analysis requires an OTU table with multiple samples")
    if(sum(abund) != "1") {abund <- tss(abund)}
    pi <- abund[abund!=0] 
    div <- sum(pi^qvalue)^(1/(1-qvalue))
    return(div) 
    }

    #If input data is an OTU table
    if(is.null(dim(abund)) == FALSE){
    if(type == "abundance"){
     if(sum(colSums(abund)) != ncol(abund)) {abund <- tss(abund)}
     if(is.null(colnames(abund)) == TRUE){colnames(abund) <- seq(1:ncol(abund))}
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
    }else if(type == "incidence"){
    #Incidence-based  
    pi <- tss(rowSums(abund != 0))
    pi <- pi[pi!=0] 
    div <- sum(pi^qvalue)^(1/(1-qvalue))
    return(div) 
    }else{
    stop("The provided type of analysis is incorrect")
    }
   
    }
#Phylodiversity
}else{
    if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")  
    if(ape::is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")  

    #If input data is a vector
    if(is.null(dim(abund)) == TRUE){
    if(type == "incidence") stop("Incidence-based analysis requires an OTU table with multiple samples")
    if(identical(sort(names(abund)),sort(tree$tip.label)) == FALSE) stop("OTU names in the vector and tree do not match")  
    if(sum(abund) != "1") {abund <- tss(abund)}
    Li <- tree$edge.length #Get branch lengths
    ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node)) #Sum relative abundances per lineage
    ai <- unlist(lapply(ltips, function(TipVector) sum(abund[TipVector]))) #Sum relative abundances per lineage
    T <- sum(Li * ai) #Get total tree depth
    Li <- Li[ai != 0] #Remove zeros
    ai <- ai[ai != 0] #Remove zeros
    phylodiv <- sum(Li/T * ai^qvalue)^(1/(1-qvalue)) #Compute phylodiversity
    return(phylodiv)
    }

    #If input data is an OTU table
    if(is.null(dim(abund)) == FALSE){    
    if(identical(sort(rownames(abund)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")  
    if(type == "abundance"){
      #Abundance-based  
      if(sum(colSums(abund)) != ncol(abund)) {abund <- tss(abund)}
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
    }else if(type == "incidence"){
    #Incidence-based  
    pi <- tss(rowSums(abund != 0))
    Li <- tree$edge.length #Get branch lengths
    ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node)) #Sum relative abundances per lineage
    ai <- unlist(lapply(ltips, function(TipVector) sum(pi[TipVector]))) #Sum relative abundances per lineage
    T <- sum(Li * ai) #Get total tree depth
    Li <- Li[ai != 0] #Remove zeros
    ai <- ai[ai != 0] #Remove zeros
    phylodiv <- sum(Li/T * ai^qvalue)^(1/(1-qvalue)) #Compute phylodiversity
    return(phylodiv)
    }else{
    stop("The provided type of analysis is incorrect")
    }
    }
}
}
