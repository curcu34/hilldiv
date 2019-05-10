\name{div.test}
\alias{div.test}
\title{Diversity comparison test between groups of samples}
\usage{div.test(otutable,qvalue,hierarchy,tree)}
\arguments{
  \item{otutable}{A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs.}
  \item{qvalue}{A positive number (>=0). It can be an integer or contain decimals.}
  \item{hierarchy}{2-column matrix | The first column lists the sample names while the second lists the groups. If provided, group profiles are plotted instead of individual profiles.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
}
\details{Performs statistical hypothesis testing for comparing the neutral or phylogenetic (if tree is provided) diversities of a given order across groups of samples. The function first runs Shapiro and Barlett tests to assess normality and homogeneity of the data, and depending on the results performs either a Student’s t-test or a Wilcoxon Rank Sum Test if there are 2 groups, and either an ANOVA or a Kruskal-Wallis test if there are multiple groups.}
\examples{
div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
div.test(otu.table,qvalue=1,hierarchy=hierarchy.table,tree=tree)
div.test(otu.table,2,hierarchy.table,tree)
}
