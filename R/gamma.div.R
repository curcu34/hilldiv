#' Gamma diversity computation (based on Hill numbers)
#' @title Gamma diversity
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords gamma partitioning hill
#' @description Compute gamma diversity of a system from a matrix (OTU table) containing multiple samples. If a tree is provided, the computed gamma diversity accounts for the phylogenetic relations across OTUs.
#' @param otutable An OTU table (matrix/data.frame) indicating the absolute or relative OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs.
#' @param qvalue A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals.
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match.
#' @param weight A vector indicating the relative weight of each sample. The order needs to be identical to the order of the samples in the OTU table. The values need to sum up to 1. If empty, all samples are weighed the same.
#' @return A gamma diversity value.
#' @usage gamma.div(otutable,qvalue,tree,weight)
#' @seealso \code{\link{div.part}}, \code{\link{alpha.div}}, \code{\link{match.data}}
#' @examples
#' gamma.div(otutable=otu.table,qvalue=1,weight=weight.vector)
#' gamma.div(otutable=otu.table,qvalue=1,tree=tree)
#' gamma.div(otu.table,1,tree,weight.vector)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' Chao, A., Chiu, C.‐H., & Hsieh, T. C. (2012). Proposing a resolution to de‐ bates on diversity partitioning. Ecology, 93, 2037–2051
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439.
#' @export gamma.div

gamma.div <- function(otutable,qvalue,tree,weight){

#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999} # change q to the limit of the unity (0.99999) if q=1
if(missing(weight)) { weight= rep(1/ncol(otutable),ncol(otutable))}

if(missing(tree)){
    #Neutral
    pi <- as.data.frame(otutable[apply(otutable, 1, function(z) !all(z==0)),]) #remove OTUs without abundances (=all-zero rows)
    pi.w <- sweep(pi,2,weight,"*") #apply weights
    div <- sum(rowSums(pi.w)^qvalue)^(1/(1-qvalue)) #apply gamma diversity formula
    return(div) #print the result
}else{
    #Phylogenetic
    if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")
    if(ape::is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")
    if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
    otutable <- as.data.frame(otutable)
    wj <- weight
    N <- ncol(otutable)
    Li <- tree$edge.length
    ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
    aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(otutable[TipVector,]))), ncol = N, byrow = TRUE)
    aij.wj <- sweep(aij, 2, wj, "*")
    ai <- rowSums(aij.wj)
    T <- sum(sweep(aij.wj, 1, Li, "*"))
    L <- matrix(rep(Li, N), ncol = N)
    Li <- Li[ai != 0] #Remove zeros
    ai <- ai[ai != 0] #Remove zeros
    wm <-  matrix(rep(wj, length(Li)), ncol = N, byrow=TRUE)
    phylodiv <- (sum(Li * (ai/T)^qvalue)^(1/(1 - qvalue)))/T
    return(phylodiv)
}
}
