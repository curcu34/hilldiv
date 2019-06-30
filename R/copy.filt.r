copy.filt <- function(abund,threshold){
  if(is.null(dim(abund)) == FALSE){
  #It is an OTU table (multiple samples)
  sapply(1:ncol(abund), function(colnum){temp = abund[,colnum]
    if(threshold == round(threshold)){
    #Absolute
    rownums = which(temp < threshold)
    }else{
    #Relative
    rownums = which(temp < sum(temp)*threshold)
    }
  abund[rownums, colnum] <<- 0})
  }else{
  #It is a vector (single sample)
    if(threshold == round(threshold)){
    #Absolute
    abund[abund < threshold] <- 0
    }else{
    #Relative
    abund[abund < sum(abund)*threshold] <- 0
    }
  }
  return(abund)
}
