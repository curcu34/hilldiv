alpha.div <- function(otutable,weight,qvalue){
    
#Quality-check
if (missing(otutable)) stop("OTU table is missing.")    
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(weight)) {weight= rep(1/ncol(otutable),ncol(otutable))}
if(missing(qvalue)) stop("q value is missing")

#Function    
if (qvalue==1) {qvalue=0.99999} # change q to the limit of the unity (0.99999) if q=1
pi <- as.data.frame(otutable[apply(otutable, 1, function(z) !all(z==0)),]) #remove OTUs without abundances (=all-zero rows) 
pi.w <- sweep(pi,2,weight,"*") #apply weights 
pi.w.q <- pi.w^qvalue #apply order of diversity 
pi.w.q[!pi] <- 0 #only necessary when q=0, to turn 1s created when making the exponentials of 0 to 0  
N <- length(weight) #calculate number of samples
div <- sum(rowSums(pi.w.q))^(1/(1-qvalue))/N #apply alpha diversity formula
return(div) #print the result
}
