hill.intext <- function(otutable,qvalue,hierarchy,tree,type,output,colour,size){

if((qvalue != 0) & (qvalue != 1) & (qvalue != 2))  stop("The order of diversity (q) must to be 0, 1 or 2.")
if(missing(output)){output="diversity"}
if(missing(type)){type="incidence"}
  
#############################  
#############################
       # INCIDENCE  #
#############################
#############################
if(type == "incidence"){

#Generate lists
if(!missing(hierarchy)){
lists <- to.inext(otutable,hierarchy)
maxsize <- max(table(hierarchy[,2]))
}else{
lists <- to.inext(otutable)
maxsize <- ncol(otutable)
}

#Calculate size
if(missing(size)){size=seq(1,maxsize*3,round(maxsize*3/20))}
                   
#############################
# NEUTRAL DIVERSITY (iNEXT) #
#############################
                   
if(missing(tree)){
  
#Run iNEXT
sp.inext <- iNEXT(lists, q=qvalue, datatype="incidence_raw", size=size)

#Return iNEXT object
if (output == "report"){
return(sp.inext)
}
  
#Extract data from iNEXT object
table <- c()
for(subsystem in c(1:length(lists))){
row <- cbind(rep(names(sp.inext$iNextEst[subsystem]),nrow((sp.inext$iNextEst[[subsystem]]))),sp.inext$iNextEst[[subsystem]][,1],sp.inext$iNextEst[[subsystem]][,2],sp.inext$iNextEst[[subsystem]][,4],sp.inext$iNextEst[[subsystem]][,5],sp.inext$iNextEst[[subsystem]][,6],sp.inext$iNextEst[[subsystem]][,7],sp.inext$iNextEst[[subsystem]][,8],sp.inext$iNextEst[[subsystem]][,9])
table <- rbind(table,row)
}
melted.inext <- as.data.frame(table)
melted.inext[,2] <- as.numeric(as.character(melted.inext[,2]))
melted.inext[,4] <- as.numeric(as.character(melted.inext[,4]))
melted.inext[,5] <- as.numeric(as.character(melted.inext[,5]))
melted.inext[,6] <- as.numeric(as.character(melted.inext[,6]))
melted.inext[,7] <- as.numeric(as.character(melted.inext[,7]))
melted.inext[,8] <- as.numeric(as.character(melted.inext[,8]))
melted.inext[,9] <- as.numeric(as.character(melted.inext[,9]))
colnames(melted.inext) <- c("Subsystem","Size","Method","Diversity","Div_min","Div_max","Completeness","Com_min","Com_max")

#Plot diversity
if (output == "diversity"){

if(missing(colour)){
getPalette <- colorRampPalette(brewer.pal(length(lists), "Paired"))
colour <- getPalette(length(lists))
}
  
plot <- ggplot() +
geom_line(data = melted.inext[which(melted.inext$Method %in% c("interpolated","observed")),], aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_line(data = melted.inext[which(melted.inext$Method  %in% c("extrapolated","observed")),], aes(x = Size, y = Diversity, colour=Subsystem), linetype=2) + 
geom_point(data = melted.inext[which(melted.inext$Method == "observed"),],aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_ribbon(data = melted.inext,aes(x = Size, ymin = Div_min, ymax = Div_max, group=Subsystem, fill=Subsystem), alpha = 0.05) +
xlab("Sample size") + 
ylab("Diversity") +
scale_colour_manual(values = colour) + 
scale_fill_manual(values = colour) + 
theme_minimal()
print(plot)
}
 
#Plot completeness
if(output == "completeness"){

if(missing(colour)){
getPalette <- colorRampPalette(brewer.pal(length(lists), "Paired"))
colour <- getPalette(length(lists))
}
  
plot <- ggplot() +
geom_line(data = melted.inext[which(melted.inext$Method %in% c("interpolated","observed")),], aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_line(data = melted.inext[which(melted.inext$Method  %in% c("extrapolated","observed")),], aes(x = Size, y = Diversity, colour=Subsystem), linetype=2) + 
geom_point(data = melted.inext[which(melted.inext$Method == "observed"),],aes(x = Size, y = Completeness, colour=Subsystem)) +
geom_ribbon(data = melted.inext,aes(x = Size, ymin = Com_min, ymax = Com_max, group=Subsystem, fill=Subsystem), alpha = 0.05) +
xlab("Sample size") + 
ylab("Completeness") +
scale_colour_manual(values = colour) + 
scale_fill_manual(values = colour) + 
theme_minimal()
print(plot)
}
  
}else if (type == "abundance_alpha"){
 
####################################
# PHYLOGENETIC DIVERSITY (iNextPD) #
####################################
  
if(class(tree) != "phylog"){tree <- phylo.to.phylog(tree)}
#Run iNextPD
sp.inextpd <- iNextPD(lists, q=qvalue, datatype="incidence_raw",phy=tree,labels=rownames(otutable),size=size)

#Return iNextPD object
if(output == "report"){
return(sp.inextpd)
} 
  
#Extract data from iNextPD object
table <- c()
for(subsystem in c(1:length(lists))){
row <- cbind(rep(names(sp.inextpd$iNextPDEst[subsystem]),nrow((sp.inextpd$iNextPDEst[[subsystem]]))),sp.inextpd$iNextPDEst[[subsystem]][,1],sp.inextpd$iNextPDEst[[subsystem]][,2],sp.inextpd$iNextPDEst[[subsystem]][,4],sp.inextpd$iNextPDEst[[subsystem]][,5])
table <- rbind(table,row)
}
melted.inextpd <- as.data.frame(table)
melted.inextpd[,2] <- as.numeric(as.character(melted.inextpd[,2]))
melted.inextpd[,4] <- as.numeric(as.character(melted.inextpd[,4]))
melted.inextpd[,5] <- as.numeric(as.character(melted.inextpd[,5]))
colnames(melted.inextpd) <- c("Subsystem","Size","Method","Diversity","Completeness")

#Plot diversity
if (output == "diversity"){

if(missing(colour)){
getPalette <- colorRampPalette(brewer.pal(length(lists), "Paired"))
colour <- getPalette(length(lists))
}  
  
plot <- ggplot() +
geom_line(data = melted.inextpd[which(melted.inextpd$Method %in% c("interpolated","observed")),], aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_line(data = melted.inextpd[which(melted.inextpd$Method  %in% c("extrapolated","observed")),], aes(x = Size, y = Diversity, colour=Subsystem), linetype=2) + 
geom_point(data = melted.inextpd[which(melted.inextpd$Method == "observed"),],aes(x = Size, y = Diversity, colour=Subsystem)) +
xlab("Sample size") + 
ylab("Diversity") +
scale_colour_manual(values = colour) + 
scale_fill_manual(values = colour) + 
theme_minimal()
print(plot)
}
 
#Plot completeness
if(output == "completeness"){

if(missing(colour)){
getPalette <- colorRampPalette(brewer.pal(length(lists), "Paired"))
colour <- getPalette(length(lists))
} 
  
plot <- ggplot() +
geom_line(data = melted.inextpd[which(melted.inextpd$Method %in% c("interpolated","observed")),], aes(x = Size, y = Diversity, colour=Subsystem)) +
geom_line(data = melted.inextpd[which(melted.inextpd$Method  %in% c("extrapolated","observed")),], aes(x = Size, y = Diversity, colour=Subsystem), linetype=2) + 
geom_point(data = melted.inextpd[which(melted.inextpd$Method == "observed"),],aes(x = Size, y = Completeness, colour=Subsystem)) +
xlab("Sample size") + 
ylab("Completeness") +
scale_colour_manual(values = colour) + 
scale_fill_manual(values = colour) + 
theme_minimal()
print(plot)
}
  
}
}
  
#############################  
#############################
     # ABUNDANCE ALPHA #
#############################
#############################
if(type == "abundance"){
  
maxsize=max(colSums(otutable))  
if(missing(size)){size=seq(1,maxsize*2,round(maxsize*2/20))}

#############################
# NEUTRAL DIVERSITY (iNEXT) #
#############################
                   
if(missing(tree)){  

#Run iNext
sp.inext <- iNEXT(otutable, q=qvalue, datatype="abundance", size=size)
                  
#Return iNEXT object
if(output == "report"){
return(sp.inext)
}
  
}else{
####################################
# PHYLOGENETIC DIVERSITY (iNextPD) #
####################################

if(class(tree) != "phylog"){tree <- phylo.to.phylog(tree)}
  
#Run iNextPD
sp.inextpd <- iNextPD(lists, q=qvalue, datatype="abundance",phy=tree,labels=rownames(otutable),size=size)

#Return iNextPD object
if(output == "report"){
return(sp.inextpd)
} 
  
  
}
  
  
}else{
#############################  
#############################
     # ABUNDANCE GAMMA #
#############################
#############################

#############################
# NEUTRAL DIVERSITY (gamma.est) #
#############################
                   
if(missing(tree)){  

#Run Gamma.est
gamma.est <- gamma.est(otutable=otutable,qvalue=qvalue,hierarchy=hierarchy,steps=100,iter=30,pred.size=40,summary=TRUE,conf=FALSE)
                  
#Return Gamma.est object
if(output == "report"){
return(gamma.est)
}
}else{
  
####################################
# PHYLOGENETIC DIVERSITY (gamma.est) #
####################################
  
#Run Gamma.est
gamma.est <- gamma.est(otutable=otutable,qvalue=qvalue,hierarchy=hierarchy,tree=tree,steps=100,iter=30,pred.size=40,summary=TRUE,conf=FALSE)
#Return Gamma.est object
if(output == "report"){
return(gamma.est)
} 
}

#Plot
  
#Declare depth check function
depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(this,depth,thisdepth=thisdepth+1))))    
  }
}

if(depth(gamma.est) = 2){
#####
# Single group
#####

steps.table <- gamma.est$Steps
steps.table$Type <- c(rep("interpolated",length(steps.table$Observed[!is.na(steps.table$Observed)])-1),"observed",rep("extrapolated",length(steps.table$Observed[is.na(steps.table$Observed)])))


ggplot() +
geom_line(data = steps.table[which(steps.table$Type %in% c("interpolated","observed")),], aes(x = Step, y = Modelled)) +
geom_line(data = steps.table[which(steps.table$Type %in% c("extrapolated","observed")),], aes(x = Step, y = Modelled), linetype=2) + 
geom_point(data = steps.table[which(steps.table$Type == "observed"),],aes(x = Step, y = Modelled)) +
#geom_ribbon(data = steps.table,aes(x = Step, ymin = Com_min, ymax = Com_max, group=Subsystem, fill=Subsystem), alpha = 0.05) +
xlab("Sample size") + 
ylab("Diversity") +
scale_colour_manual(values = colour) + 
scale_fill_manual(values = colour) + 
theme_minimal()

}else if(depth(gamma.est) = 3){
#####
# Multiple groups
#####

if(missing(colour)){
getPalette <- colorRampPalette(brewer.pal(length(gamma.est), "Paired"))
colour <- getPalette(length(gamma.est))
}

steps.table <- c()
for(g in c(1:length(gamma.est))){
group.table <- gamma.est[[g]]$Steps
group.table$Subsystem <- rep(names(gamma.est)[g],nrow(group.table))
group.table$Type <- c(rep("interpolated",length(group.table$Observed[!is.na(group.table$Observed)])-1),"observed",rep("extrapolated",length(group.table$Observed[is.na(group.table$Observed)])))
steps.table <- rbind(steps.table,group.table)
}

ggplot() +
geom_line(data = steps.table[which(steps.table$Type %in% c("interpolated","observed")),], aes(x = Step, y = Modelled, colour=Subsystem)) +
geom_line(data = steps.table[which(steps.table$Type  %in% c("extrapolated","observed")),], aes(x = Step, y = Modelled, colour=Subsystem), linetype=2) + 
geom_point(data = steps.table[which(steps.table$Type == "observed"),],aes(x = Step, y = Modelled, colour=Subsystem)) +
#geom_ribbon(data = steps.table,aes(x = Step, ymin = Com_min, ymax = Com_max, group=Subsystem, fill=Subsystem), alpha = 0.05) +
xlab("Sample size") + 
ylab("Diversity") +
scale_colour_manual(values = colour) + 
scale_fill_manual(values = colour) + 
theme_minimal()

}else{
warning("The provided gamma estimation object is not correct")
}
  
  
  
  
}
}
