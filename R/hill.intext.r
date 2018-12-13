hill.intext <- function(otutable,qvalue,hierarchy,tree){

#Generate lists
if(!missing(hierarchy)){
lists <- to.inext(otutable,hierarchy)
maxsize <- max(table(hierarchy[,2]))
}else{
lists <- to.inext(otutable)
maxsize <- ncol(otutable)
}


#Run iNEXT/iNextPD
if(!missing(tree)){
sp.inext <- iNEXT(lists, q=qvalue, datatype="incidence_freq",size=maxsize*3)
}else{
if(class(tree) != "phylog"){tree.phylog <- phylo.to.phylog(tree)}
sp.inext <- iNextPD(lists, q=qvalue, datatype="incidence_freq",phy=tree.phylog)
}

table <- c()
for(subsystem in c(1:length(lists))){
row <- cbind(rep(names(sp.inext$iNextEst[subsystem]),nrow((sp.inext$iNextEst[[subsystem]]))),sp.inext$iNextEst[[subsystem]][,1],sp.inext$iNextEst[[subsystem]][,4])
table <- rbind(table,row)
}
melted.inext <- as.data.frame(table)
melted.inext[,2] <- as.numeric(as.character(melted.inext[,2]))
melted.inext[,3] <- as.numeric(as.character(melted.inext[,3]))
colnames(melted.inext) <- c("Subsystem","Size","Diversity")

#Plot
getPalette = colorRampPalette(brewer.pal(length(lists), "Paired"))
plot <- ggplot(melted.inext , aes(x = Size, y = Diversity, group=Subsystem, colour=Subsystem)) +
geom_line() + 
xlab("Sample size") + 
ylab("Diversity") +
scale_colour_manual(values = getPalette(length(lists))) + 
theme_minimal()
print(plot)

}
