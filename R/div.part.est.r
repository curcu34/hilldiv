div.part.est <- function(otutable,qvalue,hierarchy,tree){

#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if((qvalue != 0) & (qvalue != 1) & (qvalue != 2))  stop("The order of diversity (q) must to be 0, 1 or 2.")

#L2 data processing
otutable.L2 <- merge(t(otutable),hierarchy, by.x="row.names",by.y="L1")
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- aggregate(otutable.L2[,-ncol(otutable.L2)], by=list(otutable.L2[,ncol(otutable.L2)]), function(x) sum(x != 0))
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- t(sweep(otutable.L2, 1, rowSums(otutable.L2), FUN="/"))  

#Alpha diversity (incidence-based)  
  
#Gamma diversity (incidence-based)
  

  
  
}
