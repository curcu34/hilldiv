tss <- function(abund){
  if(is.null(dim(abund)) == "TRUE"){ 
  abund.norm <- abund/sum(abund)}
  if(is.null(dim(abund)) == "FALSE"){
  abund.norm <- sweep(abund, 2, colSums(abund), FUN="/")}
  return(abund.norm)
}
