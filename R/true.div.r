true.div <- function(x,q){ 
if(q==1){q=0.99999} 
pi <- x[x!=0] 
div <- sum(pi^q)^(1/(1-q))
return(div) 
}
