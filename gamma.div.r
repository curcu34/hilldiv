gamma.div <- function(x,w,q){
if (q==1) {q=0.99999} # change q to the limit of the unity (0.99999) if q=1
pi <- as.data.frame(x[apply(x, 1, function(z) !all(z==0)),]) #remove OTUs without abundances (=all-zero rows) 
pi.w <- sweep(pi,2,w,"*") #apply weights 
div <- sum(rowSums(pi.w)^q)^(1/(1-q)) #apply alpha diversity formula
return(div) #print the result
}
