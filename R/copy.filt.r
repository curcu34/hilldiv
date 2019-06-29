copy.filt <- function(abund,threshold){
  if(is.null(dim(abund)) == FALSE){
  #It is an OTU table (multiple samples)
  sapply(1:ncol(abund), function(colnum){temp = abund[,colnum]
  rownums = which(temp < sum(temp)*threshold)
  abund[rownums, colnum] <<- 0})
  }else{
  #It is a vector (single sample)
  abund[abund < sum(abund)*threshold] <- 0
  }
  return(abund)
}
