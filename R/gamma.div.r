gamma.div <- function(otutable,qvalue,weight){

#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(weight)) { weight= rep(1/ncol(otutable),ncol(otutable))}
if(missing(weight)) warning("Assuming equal weights")
  
if (qvalue==1) {qvalue=0.99999} # change q to the limit of the unity (0.99999) if q=1
pi <- as.data.frame(otutable[apply(otutable, 1, function(z) !all(z==0)),]) #remove OTUs without abundances (=all-zero rows) 
pi.w <- sweep(pi,2,weight,"*") #apply weights 
div <- sum(rowSums(pi.w)^qvalue)^(1/(1-qvalue)) #apply alpha diversity formula
return(div) #print the result
}
