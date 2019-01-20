div.profile <- function(abund,qvalues,tree,values,hierarchy,level,colour,log){
    
#Quality-check and warnings
if(missing(abund)) stop("The abundance data is missing")
if(missing(qvalues)) {qvalues= seq(from = 0, to = 5, by = (0.1))}
if(missing(values)) {values= "FALSE"}
if(missing(level)) {level= "gamma"}
if(missing(log)) {log= "FALSE"}

#If input data is a vector
if(is.null(dim(abund)) == TRUE){
    profile <- c()
    for (o in qvalues){
    if(missing(tree)){ 
        div.value <- hilldiv::hill.div(abund,o)
        }else{
        div.value <- hilldiv::hill.div(abund,o,tree)
    }
    profile <- c(profile,div.value)
    }
    profile.melted <- as.data.frame(cbind(qvalues,profile))
    colnames(profile.melted) <- c("Order","Profile")
    plot <- ggplot(profile.melted , aes(x = Order, y = Profile)) +
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
    if(missing(hierarchy)){
        for (o in qvalues){
            if(missing(tree)){ 
            div.values <- hilldiv::hill.div(abund,o)
            }else{
            div.values <- hilldiv::hill.div(abund,o,tree)
            }
        profile <- rbind(profile,div.values)
        }
        rownames(profile) <- qvalues
        profile.melted <- as.data.frame(melt(profile))
        colnames(profile.melted) <- c("Order","Sample","Value")
        profile.melted[,1] <- as.numeric(as.character(profile.melted[,1]))
        profile.melted[,3] <- as.numeric(as.character(profile.melted[,3]))
        if(log == "TRUE"){profile.melted[,3] <- log(profile.melted[,3])}
        
        #Declare colours
	if(missing(colour) || (length(colour) != ncol(abund))){
        getPalette <- colorRampPalette(brewer.pal(ncol(abund), "Paired"))
        colour <- getPalette(ncol(abund))
        }
        
        #Plot
        plot <- ggplot(profile.melted , aes(x = Order, y = Value, group=Sample, colour=Sample)) +
        geom_line() + 
        xlab("Order of diversity") + 
        ylab(if(log == "TRUE"){"Effective number of OTUs (log-transformed)" }else{"Effective number of OTUs"}) +
        scale_colour_manual(values = colour) + 
        theme_minimal()
        print(plot)
        
    }else{
    colnames(hierarchy) <- c("Sample","Group")
    groups <- sort(unique(hierarchy$Group))   
        for (g in groups){
            samples <- as.character(hierarchy[which(hierarchy$Group == g),1])
            abund.subset <- abund[,samples]
            abund.subset <- as.data.frame(abund.subset[apply(abund.subset, 1, function(z) !all(z==0)),])
            if(!missing(tree)){                                       
            missing.otus <- setdiff(tree$tip.label,rownames(abund.subset))
            tree.subset <- drop.tip(tree,missing.otus)
            }
                 for (o in qvalues){
                        if(missing(tree)){ 
                            if(level == "gamma"){div.value <- hilldiv::gamma.div(abund.subset,o)}
                            if(level == "alpha"){div.value <- hilldiv::alpha.div(abund.subset,o)}
                            if(level == "incidence"){div.value <- hilldiv::hill.div(rowSums(abund.subset != 0)/sum(rowSums(abund.subset != 0)),o)}
                        }else{
                            if(level == "gamma"){div.value <- hilldiv::gamma.div(abund.subset,o,tree.subset)}
                            if(level == "alpha"){div.value <- hilldiv::alpha.div(abund.subset,o,tree.subset)}
                            if(level == "incidence"){div.value <- hilldiv::hill.div(rowSums(abund.subset != 0)/sum(rowSums(abund.subset != 0)),o,tree.subset)}
                        }
                 profile <- rbind(profile,cbind(g,div.value))
                 }    
          }
                                                             
    profile <- as.data.frame(cbind(profile,rep(qvalues,length(groups))))
    profile[,2] <- as.numeric(as.character(profile[,2]))
    profile[,3] <- as.numeric(as.character(profile[,3]))
    colnames(profile) <- c("Group","Value","Order")
    if(log == "TRUE"){profile[,2] <- log(profile[,2])}
        
    #Declare colours
	if(missing(colour) || (length(colour) != length(groups))){
    getPalette <- colorRampPalette(brewer.pal(length(groups), "Paired"))
    colour <- getPalette(length(groups))
    }                                                             
    
    #Plot                                                         
    plot <- ggplot(profile , aes(x = Order, y = Value, group=Group, colour=Group)) +
    geom_line() + 
    xlab("Order of diversity") + 
    ylab(if((log == "TRUE") & missing(tree)){"Effective number of OTUs (log-transformed)"}else if((log == "TRUE") & !missing(tree)){"Effective number of lineages (log-transformed)"}else if((log == "FALSE") & !missing(tree)){"Effective number of lineages"}else{"Effective number of OTUs"}) +
    scale_colour_manual(values = colour + 
    theme_minimal()
    print(plot)
                                                                                                                     
    }
                                                     

  
    
    if(values == "TRUE"){
    return(profile)
    }
                                     
    
}

}
