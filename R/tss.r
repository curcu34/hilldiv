tss <- function(abund){ 
  #If input data is a vector
  if(is.null(dim(abund)) == "TRUE"){ 
  abund.norm <- abund/sum(abund)}
#If input data is an OTU table
  if(is.null(dim(abund)) == "FALSE"){
  abund.norm <- sweep(abund, 2, colSums(abund), FUN="/")}
  return(abund.norm)
}
