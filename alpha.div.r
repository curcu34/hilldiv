alpha.div <- function(x,w,q){
if (q==1) {q=0.99999}
pi.w <- sweep(x,2,w,"*")
pi.w.q <- pi.w^q
pi.w.q[!x] <- 0
N <- length(w)
div <- sum(rowSums(pi.w.q))^(1/(1-q))/N
return(div)
}
