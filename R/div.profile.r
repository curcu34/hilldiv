div.profile <- function(abund,order){
    
#Quality-check and warnings
if(missing(abund)) stop("The abundance data is missing")
if(missing(order)) { order= seq(from = 0, to = 5, by = (0.1))}

#If input data is a vector
if(is.null(dim(abund)) == TRUE){
    profile <- c()
    for (o in order){
    div.value <- true.div(abund,o)
    profile <- c(profile,div.value)
    }
    profile.melted <- as.data.frame(cbind(order,profile))
    plot <- ggplot(profile.melted , aes(x = order, y = profile)) +
           geom_line() + 
           xlab("Order of diversity") + ylab("Effective number of OTUs") +
           theme_minimal()
    print(plot)
}
  
#If input data is an OTU table
if(is.null(dim(abund)) == FALSE){
    
    if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
    if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
    
    profile <- c()
    for (o in order){
    div.values <- true.div(abund,o)
    profile <- rbind(profile,div.values)
    }    
    rownames(profile) <- order
    profile.melted <- melt(profile)
    colnames(profile.melted) <- c("Order","Sample","Value")
    
    getPalette = colorRampPalette(brewer.pal(ncol(abund), "Paired"))
    plot <- ggplot(profile.melted , aes(x = Order, y = Value, group=Sample, colour=Sample)) +
        geom_line() + 
        xlab("Order of diversity") + ylab("Effective number of OTUs") +
        scale_colour_manual(values = getPalette(ncol(abund))) + 
        theme_minimal()
    print(plot)
}

}
