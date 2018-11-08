alpha.div <- function(x,w,q){
if (q==1) {q=0.99999}
pi <- as.data.frame(x[apply(x, 1, function(z) !all(z==0)),])
pi.w <- sweep(pi,2,w,"*")
pi.w.q <- pi.w^q
pi.w.q[!pi] <- 0
N <- length(w)
div <- sum(rowSums(pi.w.q))^(1/(1-q))/N
return(div)
}
