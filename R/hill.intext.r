hill.intext <- function(otutable,qvalue,hierarchy,tree){

#qvalue can only be 0, 1 or 2
#Need to specify observed value

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
sp.inext <- iNEXT(lists, q=qvalue, datatype="incidence_freq",size=seq(1,maxsize*3,round(maxsize*3/20)))
}else{
if(class(tree) != "phylog"){tree.phylog <- phylo.to.phylog(tree)}
sp.inext <- iNextPD(lists, q=qvalue, datatype="incidence_freq",phy=tree.phylog,size=seq(1,maxsize*3,round(maxsize*3/20)))
}

#Extract data from iNEXT object
table <- c()
for(subsystem in c(1:length(lists))){
row <- cbind(rep(names(sp.inext$iNextEst[subsystem]),nrow((sp.inext$iNextEst[[subsystem]]))),sp.inext$iNextEst[[subsystem]][,1],sp.inext$iNextEst[[subsystem]][,2],sp.inext$iNextEst[[subsystem]][,4])
table <- rbind(table,row)
}
melted.inext <- as.data.frame(table)
melted.inext[,2] <- as.numeric(as.character(melted.inext[,2]))
melted.inext[,4] <- as.numeric(as.character(melted.inext[,4]))
colnames(melted.inext) <- c("Subsystem","Size","Method","Diversity")

#Plot
getPalette = colorRampPalette(brewer.pal(length(lists), "Paired"))
plot <- ggplot() +
geom_line(data = melted.inext[which(melted.inext$Method == "interpolated"),], aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_line(data = melted.inext[which(melted.inext$Method == "extrapolated"),], aes(x = Size, y = Diversity, colour=Subsystem), linetype=2) + 
geom_point(data = melted.inext[which(melted.inext$Method == "observed"),],aes(x = Size, y = Diversity, colour=Subsystem)) +
xlab("Sample size") + 
ylab("Diversity") +
scale_colour_manual(values = getPalette(length(lists))) + 
theme_minimal()
print(plot)

}
