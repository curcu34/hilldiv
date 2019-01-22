div.test.plot <- function(divtest,chart,colour){ 
if(missing(chart)){chart="boxplot"}
if((names(divtest)[1] != "data") & (names(divtest)[2] != "normality.pvalue")) stop("The input object does not seem to be a div.test output.")

#Get data table  
divtestdata <- divtest$data  
  
#Declare colours
if(missing(colour) || (length(colour) != divtest$groups)){
getPalette <- colorRampPalette(brewer.pal(divtest$groups, "Paired"))
colour <- getPalette(divtest$groups)
}
  
#Plot
if(chart == "boxplot"){
plot <- ggplot(divtestdata, aes(x=Group, y=Value, colour=Group, fill=Group)) + 
  geom_boxplot() +
  ylab("Effective number of OTUs") + xlab("Groups") +
  scale_colour_manual(values=colour) +
  scale_fill_manual(values=alpha(colour, 0.3)) +
  theme_minimal() +
  coord_flip()
print(plot)
}
if(chart == "jitter"){
plot <- ggplot(divtestdata, aes(x=Group, y=Value, colour=Group)) + 
  geom_jitter(width = 0.1) +
  ylab("Effective number of OTUs") + xlab("Groups") +
  scale_colour_manual(values=alpha(colour, 0.6)) +
  theme_minimal() +
  coord_flip()
print(plot)
}  
if(chart == "violin"){
plot <- ggplot(divtestdata, aes(x=Group, y=Value, colour=Group)) + 
  geom_violin() +
  ylab("Effective number of OTUs") + xlab("Groups") +
  scale_colour_manual(values=alpha(colour, 0.3)) +
  theme_minimal() +
  coord_flip()
print(plot)
}  
  
}
