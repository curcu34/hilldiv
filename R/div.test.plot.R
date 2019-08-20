#' Diversity test plotting
#' @title Diversity test plotting
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords hill numbers comparison chart
#' @description Plot of diversity comparison between groups of samples
#' @param divtest Object outputed by the div.test() function
#' @param chart Chart type, either 'box' for boxplot, 'jitter' for jitter plot or 'violin' for violin plot. chart="box"
#' @param colour The number of vector items (colours, e.g. '#34k235'), must equal the number of groups that are intended to plot.
#' @param posthoc If 'TRUE' pairwise p-values of the posthoc analyses will be ploted. It requires the div.test() object to contain posthoc results.
#' @param threshold Maximum p-value to show in pairwise posthoc results (usually 0.05, but could be any other number between 0 an 1). P-values above the threshold will not be showed.
#' @return Chart of (mean) diversities of contrasting groups with optional posthoc results.
#' @seealso \code{\link{div.test}}, \code{\link{hill.div}}, \code{\link{div.part}}
#' @examples
#' contrast.div.q0 <- div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
#' div.test.plot(contrast.div.q0,chart="jitter")
#' div.test.plot(contrast.div.q0,chart="violin")
#' div.test.plot(div.test.result,chart="jitter",posthoc=TRUE,threshold=0.5)
#' @export

div.test.plot <- function(divtest,chart,colour,posthoc,threshold){
if(missing(chart)){chart="box"}
if(missing(posthoc)){posthoc=FALSE}
if((names(divtest)[1] != "data") & (names(divtest)[2] != "normality.pvalue")) stop("The input object does not seem to be a div.test output.")
if(names(divtest)[7] != "posthoc.method") stop("The input div.test object does not seem to contain pairwise posthoc data. Re-run div.test() using 'posthoc=TRUE' argument.")

#Get data table
divtestdata <- divtest$data
divtestdata$Group <- as.factor(divtestdata$Group)
#Sort factor levels
divtestdata$Group <- factor(divtestdata$Group, levels = as.character(unique(divtestdata$Group)))

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

#Prepare pairwisetable from posthoc data
if(names(divtest)[7] == "posthoc.method"){
combinations <- matrix(gsub(" $","",gsub("^ ","",unlist(strsplit(as.character(divtest$posthoc.result[,1]), "-", fixed = TRUE)))),ncol=2,byrow=TRUE)
pvalue <- round(divtest$posthoc.result[,4],3)
pairwisetable <- as.data.frame(cbind(combinations,pvalue))
colnames(pairwisetable) <- c("group1","group2","p")
}
pairwisetable[,1] <- as.character(pairwisetable[,1])
pairwisetable[,2] <- as.character(pairwisetable[,2])
pairwisetable[,3] <- as.numeric(as.character(pairwisetable[,3]))

#Filter pairwisetable
if(!missing(threshold)){
pairwisetable <- pairwisetable[which(pairwisetable$p < threshold),]
}

#Set y values
sortedgroups <- unique(sort(c(pairwisetable$group1,pairwisetable$group2)))
datamax <- round(max(divtest$data[which(divtest$data[,3] %in% sortedgroups),2]))
datamin <- round(min(divtest$data[which(divtest$data[,3] %in% sortedgroups),2]))
datarange <- datamax - datamin
by <- datarange * 0.1
min <- datamax
max <- min + (by*nrow(pairwisetable))
ypos <- seq(min,max,by)[-1]
pairwisetable$ypos <- ypos

#Plot
if(chart == "box"){
plot <- ggboxplot(divtestdata, x='Group', y='Value', color='Group', fill='Group', x.text.angle = 45) +
      ylab("Effective number of OTUs") + xlab("Groups") +
      scale_colour_manual(values=scales::alpha(colour, 1)) +
      scale_fill_manual(values=scales::alpha(colour, 0.5))
  if(posthoc == TRUE){
    plot <- suppressWarnings(plot + stat_pvalue_manual(pairwisetable, label = "p", y.position = "ypos"))
  }
print(plot)
}

if(chart == "jitter"){
plot <- ggboxplot(divtestdata, x='Group', y='Value', color='Group', add = "jitter", width = 0, x.text.angle = 45) +
      ylab("Effective number of OTUs") + xlab("Groups") +
      scale_colour_manual(values=scales::alpha(colour, 0))
    if(posthoc == TRUE){
      plot <- suppressWarnings(plot + stat_pvalue_manual(pairwisetable, label = "p", y.position = "ypos"))
    }
  print(plot)
  }

if(chart == "violin"){
  plot <- ggviolin(divtestdata, x='Group', y='Value', color='Group', fill='Group', x.text.angle = 45) +
        ylab("Effective number of OTUs") + xlab("Groups") +
        scale_fill_manual(values=scales::alpha(colour, 0.1)) +
        scale_colour_manual(values=scales::alpha(colour, 1))
    if(posthoc == TRUE){
      plot <- suppressWarnings(plot + stat_pvalue_manual(pairwisetable, label = "p", y.position = "ypos"))
    }
  print(plot)
  }

}
