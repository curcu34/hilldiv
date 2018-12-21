CqN <- function(beta,qvalue,N){ 
value = ((1/beta)^(qvalue-1) - (1/N)^(qvalue-1)) / (1 - (1/N)^(qvalue-1))
return(value)
}
