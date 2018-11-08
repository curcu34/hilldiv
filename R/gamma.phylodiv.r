gamma.phylodiv <- function(x,w,t,q){
if (q==1) {q=0.99999}
x <- as.data.frame(x)
wj <- w
N <- ncol(x)
Li <- t$edge.length
ltips <- sapply(t$edge[, 2], function(node) geiger::tips(t, node))
aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(x[TipVector,]))), ncol = N, byrow = TRUE)
aij.wj <- sweep(aij, 2, wj, "*")
ai <- rowSums(aij.wj)
T <- sum(sweep(aij.wj, 1, Li, "*"))
L <- matrix(rep(Li, N), ncol = N)
wm <-  matrix(rep(wj, length(Li)), ncol = N, byrow=TRUE)
phylodiv <- sum(Li * (ai/T)^q)^(1/(1 - q))
return(phylodiv)
}
