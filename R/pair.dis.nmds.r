pair.dis.nmds <- function(distance,hierarchy){

matrix = pair.div.q2$L1_UqN
values.NMDS<-metaMDS(as.dist(matrix), k = 2, trymax = 400)
NMDS=data.frame(x=values.NMDS$point[,1],y=values.NMDS$point[,2],Sample=as.factor(sample.species[,1]),Species=as.factor(sample.species[,2]))

#NMDS plot
species.colours <- read.table("1-Files/species.colours.txt", sep="\t",header=TRUE, comment.char = "")
my_col_scheme <- as.character(species.colours[,2])
ggplot() + 
	geom_point(data = NMDS, aes(x=x, y=y, colour=Species), size = 2, alpha = 0.5) +
	scale_colour_manual(values = my_col_scheme) +
	scale_shape_manual(values=16) +
	theme(panel.background = element_rect(fill = 'white', colour = 'grey'))


}
