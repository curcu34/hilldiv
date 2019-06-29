copy.filt <- function(otutable,threshold){
  sapply(1:ncol(otutable), function(colnum){temp = otutable[,colnum]
  rownums = which(temp < sum(temp)*threshold)
  otutable[rownums, colnum] <<- 0})
  return(otutable)
}
