alpha.div <- function(otutable,qvalue,tree,weight){
    
#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(weight)) { weight= rep(1/ncol(otutable),ncol(otutable))}
if(missing(weight)) warning("Assuming equal weights")

#Neutral  
if(missing(tree)){ 
        if (qvalue==1) {qvalue=0.99999} # change q to the limit of the unity (0.99999) if q=1
        pi <- as.data.frame(otutable[apply(otutable, 1, function(z) !all(z==0)),]) #remove OTUs without abundances (=all-zero rows) 
        pi.w <- sweep(pi,2,weight,"*") #apply weights 
        pi.w.q <- pi.w^qvalue #apply order of diversity 
        pi.w.q[!pi] <- 0 #only necessary when q=0, to turn 1s created when making the exponentials of 0 to 0  
        N <- length(weight) #calculate number of samples
        div <- sum(rowSums(pi.w.q))^(1/(1-qvalue))/N #apply alpha diversity formula
        return(div) #print the result
}else{
#Non-neutral  
        if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
        if(ape::is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")  
        if (qvalue==1) {qvalue=0.99999}
        otutable <- as.data.frame(otutable)
        wj <- weight
        N <- ncol(otutable)
        Li <- tree$edge.length
        ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
        aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(otutable[TipVector,]))), ncol = N, byrow = TRUE)
        aij.wj <- sweep(aij, 2, wj, "*")
        T <- sum(sweep(aij.wj, 1, Li, "*"))
        L <- matrix(rep(Li, N), ncol = N)
        wm <-  matrix(rep(wj, length(Li)), ncol = N,byrow=TRUE)
        i <-  which(aij > 0)
        phylodiv <- sum(L[i] * (aij[i]*wm[i]/T)^qvalue)^(1/(1 - qvalue))/(N*T)
        return(phylodiv)
}
                                   
}
