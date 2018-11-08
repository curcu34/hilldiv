true.div <- function(vector,qvalue){ 
  
#Quality-check and warnings
if(missing(vector)) stop("Vector is missing")
if(sum(vector) != 1) stop("The vector does not sum up to 1")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
  
#Function
if(qvalue==1){qvalue=0.99999} 
pi <- vector[vector!=0] 
div <- sum(pi^qvalue)^(1/(1-qvalue))
return(div) 
}
