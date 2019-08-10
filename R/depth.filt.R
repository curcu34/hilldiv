depth.filt <- function(otutable,threshold){
if(colSums(otutable)[1] == 1)stop("OTU table contains relative abundances. This function is only meaningful for absolute abundances.")
otutable.filt <- otutable[,colSums(otutable)>threshold]
return(otutable.filt)
}
