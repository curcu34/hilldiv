true.phylodiv <- function(x,q,t){
if(q==1){q=0.99999}

#Get branch lengths
Li <- t$edge.length

#Sum relative abundances per lineage
ltips <- sapply(t$edge[, 2], function(node) geiger::tips(t, node))
ai <- unlist(lapply(ltips, function(TipVector) sum(x[TipVector])))

#Get total tree depth
T <- sum(Li * ai)

#Eliminate zeros
Li <- Li[ai != 0]
ai <- ai[ai != 0]

#Compute phylodiversity
phylodiv <- sum(Li/T * ai^q)^(1/(1-q))

return(phylodiv)
}
