div.profile <- function(abund,tree,order,values,hierarchy){
    
#Quality-check and warnings
if(missing(abund)) stop("The abundance data is missing")
if(missing(order)) { order= seq(from = 0, to = 5, by = (0.1))}
if(missing(values)) { values= "FALSE"}

#If input data is a vector
if(is.null(dim(abund)) == TRUE){
    profile <- c()
    for (o in order){
    if(missing(tree)){ 
        div.value <- hilldiv::hill.div(abund,o)
        }else{
        div.value <- hilldiv::hill.div(abund,o,tree)
    }
    profile <- c(profile,div.value)
    }
    profile.melted <- as.data.frame(cbind(order,profile))
    plot <- ggplot(profile.melted , aes(x = order, y = profile)) +
           geom_line() + 
           xlab("Order of diversity") + ylab("Effective number of OTUs") +
           theme_minimal()
    print(plot)
    if(values == "TRUE"){
    return(profile)
    }
}
  
#If input data is an OTU table
if(is.null(dim(abund)) == FALSE){
    
    if(dim(abund)[1] < 2) stop("The OTU table only less than 2 OTUs")
    if(dim(abund)[2] < 2) stop("The OTU table contains less than 2 samples")
        
    profile <- c()
    for (o in order){
        if(missing(tree)){ 
        div.values <- hilldiv::hill.div(abund,o)
        }else{
        div.values <- hilldiv::hill.div(abund,o,tree)
    }
    profile <- rbind(profile,div.values)
    }    
    rownames(profile) <- order
    profile.melted <- melt(profile)
    colnames(profile.melted) <- c("Order","Sample","Value")
    getPalette = colorRampPalette(brewer.pal(ncol(abund), "Paired"))

    if(hasArg(hierarchy) == TRUE){
    colnames(hierarchy) <- c("Sample","Group")
    profile.melted2 <- merge(profile.melted,hierarchy,by="Sample")
    profile.melted_mean <- aggregate(profile.melted2[,"Value"],by=list(profile.melted2[,"Group"],FUN=mean)
    profile.melted_sterr <- aggregate(profile.melted2[,"Value"],by=list(profile.melted2[,"Group"],FUN=sd)#change for standard error
    #add plot code with mean+stderr
    }else{    
    plot <- ggplot(profile.melted , aes(x = Order, y = Value, group=Sample, colour=Sample)) +
        geom_line() + 
        xlab("Order of diversity") + ylab("Effective number of OTUs") +
        scale_colour_manual(values = getPalette(ncol(abund))) + 
        theme_minimal()
    }
    print(plot)
    if(values == "TRUE"){
    return(profile)
    }
                                     
    
}

}
