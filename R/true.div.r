true.div <- function(abund,qvalue){ 
  
#Quality-check and warnings
if(missing(abund)) stop("Abundance data is missing")
if(missing(qvalue)) stop("q value is missing")
if(qvalue==1){qvalue=0.99999} 
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")

#If input data is a vector
if(is.null(dim(abund)) == TRUE){
if(sum(abund) != 1) stop("The abundance data does not sum up to 1")
pi <- abund[abund!=0] 
div <- sum(pi^qvalue)^(1/(1-qvalue))
return(div) 
}

#If input data is an OTU table
if(is.null(dim(abund)) == FALSE){
samples <- colnames(abund)

sample.vector <- c()
for (s in samples){
vector <- abund[,s]
pi <- vector[vector!=0] 
div <- sum(pi^qvalue)^(1/(1-qvalue))
names(div) <- s
sample.vector <- c(sample.vector,div)
}
return(sample.vector) 
}
  
}
