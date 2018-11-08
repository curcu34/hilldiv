true.phylodiv <- function(vector,qvalue,tree){
  
#Quality-check and warnings
if(missing(vector)) stop("Vector is missing")
if(sum(vector) != 1) stop("The vector does not sum up to 1")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(tree)) stop("Tree is missing")
if(identical(sort(names(vector)),sort(tree$tip.label)) == FALSE) stop("OTU names in the vector and tree do not match")

#Function
if(qvalue==1){qvalue=0.99999}
Li <- tree$edge.length #Get branch lengths
ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node)) #Sum relative abundances per lineage
ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector]))) #Sum relative abundances per lineage
T <- sum(Li * ai) #Get total tree depth
Li <- Li[ai != 0] #Remove zeros
ai <- ai[ai != 0] #Remove zeros
phylodiv <- sum(Li/T * ai^qvalue)^(1/(1-qvalue)) #Compute phylodiversity
return(phylodiv) #Return value
}
