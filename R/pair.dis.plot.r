pair.dis.plot <- function(distance,hierarchy,type,level,colour){

if(missing(type)){type = "NMDS"}	
if(missing(level)){level = "1"}	
	
values.NMDS<-metaMDS(as.dist(distance), k = 2, trymax = 400)
if(level == 1){
NMDS=data.frame(x=values.NMDS$point[,1],y=values.NMDS$point[,2],Sample=as.factor(hierarchy[,1]),Group=as.factor(hierarchy[,2]))
}
if(level == 2){	
NMDS=data.frame(x=values.NMDS$point[,1],y=values.NMDS$point[,2],Group=as.factor(unique(hierarchy[,2])))
}
	
#Declare colours
	
if(missing(colour) || (length(colour) != length(unique(hierarchy[,2])))){
getPalette <- colorRampPalette(brewer.pal(length(unique(hierarchy[,2])), "Paired"))
colour <- getPalette(length(unique(hierarchy[,2])))
}

if(type == "NMDS"){
#NMDS plot
nmds.plot <- ggplot() + 
	geom_point(data = NMDS, aes(x=x, y=y, colour=Group), size = 2, alpha = 0.5) +
	scale_colour_manual(values = colour) +
	scale_shape_manual(values=16) +
	theme(panel.background = element_rect(fill = 'white', colour = 'grey'))
print(nmds.plot)
}
if(type == "heatmap"){
#Heatmap plot
}
}
