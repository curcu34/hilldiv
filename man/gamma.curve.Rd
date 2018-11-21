\name{gamma.curve}
\alias{gamma.curve}
\title{Gamma curve}
\usage{gamma.curve(otutable,qvalue,tree,hierarchy,steps,iter,summary)}
\arguments{
  \item{otutable}{A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs. The relative abundances of each sample (columns) must sum to 1.}
  \item{qvalue}{A positive number (>0). It can be an integer or contain decimals.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the matrix rows. Note that including tree information considerably slows down the generation of the gamma curve, hence lowering the step percentage and the number of iterations is advisable.}
  \item{hierarchy}{2-column matrix | The first column lists the sample names while the second lists the groups. If provided, group gamma curves are plotted instead of gamma curves.}
  \item{steps}{ Percentage | Percentage of samples to be used to draw the curve. '100' means rarefaction values will be generated for all sample sizes from 0 to actual sample size, e.g. 1,2,3... '50' means rarefaction values will be generated every two steps (50% of possible sample sizes) from 0 to actual sample size, e.g. 1,3,5... steps=100}
  \item{iter}{ Integer | The number of iterations to create the gamma curve. The higher the number, the lower is the standard error and the smoother is the curve. iter=50}
  \item{summary}{ Logical (TRUE/FALSE) | If 'TRUE' the summary data used to generate the gamma curve is also printed. }
\details{The gamma curves show the cummulative diversity across samples.}
\examples{
gamma.curve(otutable=otu.table,qvalue=1,tree=otu.tree,hierarchy=hierarchy.table,steps=25,iter=10,summary="TRUE")
gamma.curve(otutable,1,otu.tree,hierarchy.table)
gamma.curve(otu.table,1)
}
