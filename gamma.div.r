gamma.div <- function(x,w,q){
if (q==1) {q=0.99999}
pi.w <- sweep(x,2,w,"*")
div <- sum(rowSums(pi.w)^q)^(1/(1-q))
return(div)
}
