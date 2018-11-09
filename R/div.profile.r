div.profile <- function(abund,order){
    
#Quality-check and warnings
if(missing(abund)) stop("The abundance data is missing")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(order)) { order= seq(from = 0, to = 5, by = (0.1))}

#If input data is a vector
if(is.null(dim(abund)) == TRUE){
    profile <- c()
    for (o in order){
    div.value <- true.div(abund,o)
    profile <- c(profile,div.value)
    }
    profile.melted <- as.data.frame(cbind(order,profile))
    ggplot(profile.melted , aes(x = order, y = profile)) +
           geom_line() + 
           xlab("Order of diversity") + ylab("Effective number of OTUs") +
           theme_minimal()
}
  
#If input data is an OTU table
if(is.null(dim(abund)) == FALSE){
    
    if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
    if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
    
    profile <- c()
    for (o in order){
    div.value <- true.div(abund,o)
    profile <- c(profile,div.value)
    }    
    
}

}
