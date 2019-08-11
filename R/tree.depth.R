tree.depth <- function(tree,abund){
if(is.null(dim(abund)) == TRUE){vector <- tss(abund)}  
if(is.null(dim(abund)) == FALSE){vector <- tss(abund[,1])}  
vector <- tss(vector)
Li <- tree$edge.length
ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
T <- sum(Li * ai)
return(T)
}
