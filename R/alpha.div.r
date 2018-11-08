alpha.div <- function(x,w,q){
#Quality-check
try(if(exists("x") == FALSE) stop("OTU table does not exist"))
try(if(dim(x)[1] < 2) stop("The OTU table only less than 2 OTUs"))
try(if(dim(x)[2] < 2) stop("The OTU table contains less than 2 samples"))
    
#Function    
if (q==1) {q=0.99999} # change q to the limit of the unity (0.99999) if q=1
pi <- as.data.frame(x[apply(x, 1, function(z) !all(z==0)),]) #remove OTUs without abundances (=all-zero rows) 
pi.w <- sweep(pi,2,w,"*") #apply weights 
pi.w.q <- pi.w^q #apply order of diversity 
pi.w.q[!pi] <- 0 #only necessary when q=0, to turn 1s created when making the exponentials of 0 to 0  
N <- length(w) #calculate number of samples
div <- sum(rowSums(pi.w.q))^(1/(1-q))/N #apply alpha diversity formula
return(div) #print the result
}
