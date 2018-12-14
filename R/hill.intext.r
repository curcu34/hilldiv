hill.intext <- function(otutable,qvalue,hierarchy,tree,output){

if((qvalue =! 0) | (qvalue =! 1) | qvalue =! 2))  stop("The order of diversity (q) must to be 0, 1 or 2.")
if(!missing(output)){output="diversity")
 
#Generate lists
if(!missing(hierarchy)){
lists <- to.inext(otutable,hierarchy)
maxsize <- max(table(hierarchy[,2]))
}else{
lists <- to.inext(otutable)
maxsize <- ncol(otutable)
}

#NEUTRAL DIVERSITY
if(!missing(tree)){
  
#Run iNEXT
sp.inext <- iNEXT(lists, q=qvalue, datatype="incidence_raw",size=seq(1,maxsize*3,round(maxsize*3/20)))

if (output="report"){
return(report)
}
  
#Extract data from iNEXT object
table <- c()
for(subsystem in c(1:length(lists))){
row <- cbind(rep(names(sp.inext$iNextEst[subsystem]),nrow((sp.inext$iNextEst[[subsystem]]))),sp.inext$iNextEst[[subsystem]][,1],sp.inext$iNextEst[[subsystem]][,2],sp.inext$iNextEst[[subsystem]][,4],sp.inext$iNextEst[[subsystem]][,5],sp.inext$iNextEst[[subsystem]][,6])
table <- rbind(table,row)
}
melted.inext <- as.data.frame(table)
melted.inext[,2] <- as.numeric(as.character(melted.inext[,2]))
melted.inext[,4] <- as.numeric(as.character(melted.inext[,4]))
melted.inext[,5] <- as.numeric(as.character(melted.inext[,5]))
melted.inext[,6] <- as.numeric(as.character(melted.inext[,6]))
colnames(melted.inext) <- c("Subsystem","Size","Method","Diversity","Min","Max")

#Plot diversity
if (output="diversity"){
getPalette = colorRampPalette(brewer.pal(length(lists), "Paired"))
plot <- ggplot() +
geom_line(data = melted.inext[which(melted.inext$Method == "interpolated"),], aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_line(data = melted.inext[which(melted.inext$Method == "extrapolated"),], aes(x = Size, y = Diversity, colour=Subsystem), linetype=2) + 
geom_point(data = melted.inext[which(melted.inext$Method == "observed"),],aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_ribbon(data = melted.inext,aes(x = Size, ymin = Min, ymax = Max, group=Subsystem, fill=Subsystem), alpha = 0.05) +
xlab("Sample size") + 
ylab("Diversity") +
scale_colour_manual(values = getPalette(length(lists))) + 
scale_fill_manual(values = getPalette(length(lists))) + 
theme_minimal()
print(plot)
}
  
}else{
#PHYLOGENETIC DIVERSITY
  
if(class(tree) != "phylog"){tree.phylog <- phylo.to.phylog(tree)}
#Run iNextPD
sp.inextpd <- iNextPD(lists, q=qvalue, datatype="incidence_raw",phy=tree.phylog,size=seq(1,maxsize*3,round(maxsize*3/20)))
  
#Extract data from iNextPD object
table <- c()
for(subsystem in c(1:length(lists))){
row <- cbind(rep(names(sp.inextpd$iNextPDEst[subsystem]),nrow((sp.inextpd$iNextPDEst[[subsystem]]))),sp.inextpd$iNextPDEst[[subsystem]][,1],sp.inextpd$iNextPDEst[[subsystem]][,2],sp.inextpd$iNextPDEst[[subsystem]][,4])
table <- rbind(table,row)
}
melted.inextpd <- as.data.frame(table)
melted.inextpd[,2] <- as.numeric(as.character(melted.inextpd[,2]))
melted.inextpd[,4] <- as.numeric(as.character(melted.inextpd[,4]))
colnames(melted.inextpd) <- c("Subsystem","Size","Method","Diversity")

#Plot diversity
if (output="diversity"){
getPalette = colorRampPalette(brewer.pal(length(lists), "Paired"))
plot <- ggplot() +
geom_line(data = melted.inextpd[which(melted.inextpd$Method == "interpolated"),], aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_line(data = melted.inextpd[which(melted.inextpd$Method == "extrapolated"),], aes(x = Size, y = Diversity, colour=Subsystem), linetype=2) + 
geom_point(data = melted.inextpd[which(melted.inextpd$Method == "observed"),],aes(x = Size, y = Diversity, colour=Subsystem)) +
xlab("Sample size") + 
ylab("Diversity") +
scale_colour_manual(values = getPalette(length(lists))) + 
scale_fill_manual(values = getPalette(length(lists))) + 
theme_minimal()
print(plot)
}
  
  
}

}
