tree.depth <- function(tree,vector){
vector <- tss(vector)
Li <- tree$edge.length
ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
T <- sum(Li * ai)
return(T)
}
