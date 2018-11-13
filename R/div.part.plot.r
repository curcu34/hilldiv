div.part.plot <- function(otutable,qvalue,hierarchy,tree) {
  
#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}
if(missing(hierarchy)) warning("Assuming a two-level hierarchy: 1) sample, 2) total dataset")

#Function for 2-level hierarchy
if(missing(hierarchy)){
if(missing(tree)){
  L1_div <- hilldiv::alpha.div(otutable,qvalue)
  L2_div <- hilldiv::gamma.div(otutable,qvalue)
  }else{
  L1_div <- hilldiv::alpha.div(otutable,qvalue,tree)
  L2_div <- hilldiv::gamma.div(otutable,qvalue,tree)
}
beta <- L2_div/L1_div
part.table <- as.data.frame(rbind(cbind("Sample",L1_div),cbind("Overall",L2_div)))
part.table[,1] <- factor(part.table[,1], levels=c("Overall","Sample"))
hier.n <- 2
}
  
#Function for 3-level hierarchy
if(ncol(hierarchy) == 2){
hierarchy[,1] <- as.character(hierarchy[,1])
hierarchy[,2] <- as.character(hierarchy[,2])
if(identical(sort(colnames(otutable)),sort(hierarchy[,1])) == FALSE) stop("OTU names in the OTU table and the hierarchy table do not match")
hierarchy.names <- colnames(hierarchy)
colnames(hierarchy) <- c("L1","L2")
if(missing(tree)){  
  L1_div <- hilldiv::alpha.div(otutable,qvalue)
  L3_div <- hilldiv::gamma.div(otutable,qvalue)
  }else{
  L1_div <- hilldiv::alpha.div(otutable,qvalue,tree)
  L3_div <- hilldiv::gamma.div(otutable,qvalue,tree)
}
otutable.L2 <- merge(t(otutable),hierarchy, by.x="row.names",by.y="L1")
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- aggregate(otutable.L2[,-ncol(otutable.L2)], by=list(otutable.L2[,ncol(otutable.L2)]), FUN=sum)
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- t(sweep(otutable.L2, 1, rowSums(otutable.L2), FUN="/"))
if(missing(tree)){ 
  L2_div <- hilldiv::alpha.div(otutable.L2,qvalue)
  }else{
  L2_div <- hilldiv::alpha.div(otutable.L2,qvalue,tree)
}
beta_1 <- L2_div/L1_div
beta_2 <- L3_div/L2_div
part.table <- as.data.frame(rbind(cbind(hierarchy.names[1],L1_div),cbind(hierarchy.names[2],L2_div),cbind("Overall",L3_div)))
part.table[,1] <- factor(part.table[,1], levels=c("Overall",hierarchy.names[2],hierarchy.names[1]))
hier.n <- 3
}
  
#Function for 4-level hierarchy
if(ncol(hierarchy) == 3){


} 
  
#Error if more hierarchical levels
if(ncol(hierarchy) > 3){
stop("The maximum number of hierarchical levels allowed is 4")
} 

#Prepare table and plot
part.table[,2] <- as.numeric(as.character(part.table[,2]))
colnames(part.table) <- c("Level","Value")
getPalette = colorRampPalette(brewer.pal(hier.n, "Paired"))
if(missing(tree)){  
  xlabel <- "Effective number of OTUs"
  }else{
  xlabel <- "Effective number of lineages"
}
plot <- ggplot(part.table, aes(y = Value, x = Level, fill=Level)) + 
	geom_bar(stat = "identity") +
	ylab(xlabel) + xlab("Hierarchical levels") +
	scale_fill_manual( values = getPalette(hier.n) )+
	theme_minimal() +
	coord_flip()
print(plot)
}
