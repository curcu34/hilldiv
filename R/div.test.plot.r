div.test.plot <- function(divtest,chart,colour,stat,comb,symbol,flip){
if(missing(chart)){chart="box"}
if(missing(stat)){stat=FALSE}
if(missing(symbol)){symbol=FALSE}
if(missing(flip)){flip=FALSE}
if((names(divtest)[1] != "data") & (names(divtest)[2] != "normality.pvalue")) stop("The input object does not seem to be a div.test output.")

#Get data table
divtestdata <- divtest$data

#Declare colours
if(missing(colour) || (length(colour) < divtest$groups)){
getPalette <- colorRampPalette(brewer.pal(divtest$groups, "Paired"))
colour <- getPalette(divtest$groups)
}

#Prepare tests for ggpubr
if(divtest$method == "Kruskal-Wallis Test"){
g.method="kruskal.test"
pw.method="wilcox.test"
}
if(divtest$method == "ANOVA"){
g.method="anova"
pw.method="t.test"
}
if(divtest$method == "Student's t-Test"){
g.method="t.test"
}
if(divtest$method == "Wilcoxon Rank Sum Test"){
g.method="wilcox.test"
}

#Create combination list for ggpubr
if(missing(comb)){
comb <- list()
comb.each <- combn(unique(div.test.result$data$Group),m=2)
for (i in c(1:ncol(comb.each))){
comb[[i]] <- as.character(comb.each[,i])
}
}

#Plot
if(chart == "box"){
plot <- ggplot(divtestdata, aes(x=Group, y=Value, colour=Group, fill=Group)) +
  geom_boxplot() +
  ylab("Effective number of OTUs") + xlab("Groups") +
  scale_colour_manual(values=colour) +
  scale_fill_manual(values=scales::alpha(colour, 0.3)) +
  theme_minimal()
  if(stat == TRUE){
    if(symbol == TRUE){
    plot <- plot + stat_compare_means(comparisons = comb, method = pw.method, label = "p.signif")
    }else{
    plot <- plot + stat_compare_means(comparisons = comb, method = pw.method)
    }
  }
  if(flip == TRUE){
  plot <- plot + coord_flip()
  }
print(plot)
}

if(chart == "jitter"){
plot <- ggplot(divtestdata, aes(x=Group, y=Value, colour=Group)) +
  geom_jitter(width = 0.1) +
  ylab("Effective number of OTUs") + xlab("Groups") +
  scale_colour_manual(values=scales::alpha(colour, 0.6)) +
  theme_minimal()
  if(stat == TRUE){
    if(symbol == TRUE){
    plot <- plot + stat_compare_means(comparisons = comb, method = pw.method, label = "p.signif")
    }else{
    plot <- plot + stat_compare_means(comparisons = comb, method = pw.method)
    }
  }
  if(flip == TRUE){
  plot <- plot + coord_flip()
  }
  print(plot)
  }

if(chart == "violin"){
plot <- ggplot(divtestdata, aes(x=Group, y=Value, colour=Group, fill=Group)) +
  geom_violin() +
  ylab("Effective number of OTUs") + xlab("Groups") +
  scale_colour_manual(values=colour) +
  scale_fill_manual(values=scales::alpha(colour, 0.3)) +
  theme_minimal()
  if(stat == TRUE){
    if(symbol == TRUE){
    plot <- plot + stat_compare_means(comparisons = comb, method = pw.method, label = "p.signif")
    }else{
    plot <- plot + stat_compare_means(comparisons = comb, method = pw.method)
    }
  }
  if(flip == TRUE){
  plot <- plot + coord_flip()
  }
  print(plot)
  }

}
