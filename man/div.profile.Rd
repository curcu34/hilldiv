\name{div.profile}
\alias{div.profile}
\title{Diversity profile}
\usage{div.profile(abund,tree,order,value)}
\arguments{
  \item{abund}{A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs. The relative abundances of each sample (columns) must sum to 1.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{order}{A vector of sequential orders of diversity (usually from 0 to 5). order=seq(from = 0, to = 5, by = (0.1))}
  \item{value}{Logical (TRUE/FALSE). If 'TRUE' the Hill numbers for each order of diversity are printed on screen. value="FALSE"}
}
\details{The diversity profile enables visualizing the diversity distribution of a sample or a group of samples at a glance.}
\examples{
div.profile(abund=otu.table[,1],tree=tree,order=seq(from = 0, to = 5, by = (0.1)),value="TRUE")
div.profile(vector)
div.profile(otu.table,tree)
}
