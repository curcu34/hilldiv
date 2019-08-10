UqN <- function(beta,qvalue,N){ 
value = ((1/beta)^(1-qvalue) - (1/N)^(1-qvalue)) / (1 - (1/N)^(1-qvalue))
return(value)
}
