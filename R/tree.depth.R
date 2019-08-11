tree.depth <- function(tree,abund){
if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")  
if(is.null(dim(abund)) == TRUE){
  vector <- abund
  }  
if(is.null(dim(abund)) == FALSE){
  vector <- abund[,1]
  names(vector) <- rownames(abund)
  } 
if(identical(sort(names(vector)),sort(tree$tip.label)) == FALSE) stop("OTU names in the vector and tree do not match")  
vector <- tss(vector)
Li <- tree$edge.length
ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
T <- sum(Li * ai)
return(T)
}
