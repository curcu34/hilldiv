div.comp.plot <- function(otutable,qvalue,hierarchy,tree,chart){ 
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(hierarchy)) stop("Hierarchy table is necessary to contrast groups of samples")
if(missing(chart)) {chart="boxplot"}

if(missing(tree)){
div.values <- hill.div(otutable,qvalue)
}else{
if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")  
if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
div.values <- hill.div(otutable,qvalue,tree)
}
div.values.groups <- merge(t(t(div.values)),hierarchy,by.x="row.names",by.y="Sample")
colnames(div.values.groups) <- c("Sample","Value","Group")

#Plot
if(chart == "boxplot"){
getPalette = colorRampPalette(brewer.pal(length(unique(div.values.groups$Group)), "Paired"))
plot <- ggplot(div.values.groups, aes(x=Group, y=Value, colour=Group)) + 
  geom_boxplot() +
  xlab("Effective number of OTUs") + ylab("Groups") +
  theme_minimal() +
  coord_flip()
print(plot)
}
if(chart == "jitter"){
getPalette = colorRampPalette(brewer.pal(length(unique(div.values.groups$Group)), "Paired"))
plot <- ggplot(div.values.groups, aes(x=Group, y=Value, colour=Group)) + 
  geom_jitter(width = 0.1) +
  xlab("Effective number of OTUs") + ylab("Groups") +
  theme_minimal() +
  coord_flip()
print(plot)
}  
  
}
