#' Diversity test plotting
#' @title Diversity test plotting
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords hill numbers comparison chart
#' @description Plot of diversity comparison between groups of samples
#' @param divtest Object outputed by the div.test() function
#' @param chart Chart type, either 'box' for boxplot, 'jitter' for jitter plot or 'violin' for violin plot. chart="box"
#' @param colour The number of vector items (colours, e.g. '#34k235'), must equal the number of groups that are intended to plot.
#' @param stat If 'TRUE' pairwise mean comparison significance values will be ploted.
#' @param comb List of pairwise combinations. If missing all combinations will be plotted. e.g. list(c("Myotis myotis","Myotis capaccinii"),c("Myotis myotis","Myotis daubentonii")).
#' @param symbol If 'TRUE' symbols rather than p-values will be shown. ns: p > 0.05; *: p <= 0.05; **: p <= 0.01; ***: p <= 0.001; ****: p <= 0.0001.
#' @param flip If 'TRUE' the chart will be flipped 90 degrees.
#' @return Chart of (mean) diversities of contrasting groups.
#' @seealso \code{\link{div.test}}, \code{\link{hill.div}}, \code{\link{div.part}}
#' @examples
#' contrast.div.q0 <- div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
#' div.test.plot(contrast.div.q0,chart="jitter")
#' div.test.plot(contrast.div.q0,chart="violin")
#' div.test.plot(div.test.result,chart="jitter",stat=TRUE,flip=TRUE)
#' div.test.plot(div.test.result,stat=TRUE,comb=list(c("Myotis myotis","Myotis capaccinii")),symbol=TRUE)
#' @export

div.test.plot <- function(divtest,chart,colour,stat,comb,symbol,flip){
if(missing(chart)){chart="box"}
if(missing(stat)){stat=FALSE}
if(missing(symbol)){symbol=FALSE}
if(missing(flip)){flip=FALSE}
if((names(divtest)[1] != "data") & (names(divtest)[2] != "normality.pvalue")) stop("The input object does not seem to be a div.test output.")

#Get data table
divtestdata <- divtest$data
divtestdata$Group <- as.factor(divtestdata$Group)

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
comb.each <- combn(unique(divtestdata$Group),m=2)
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
