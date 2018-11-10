
\name{true.div}
\alias{true.div}
\title{True diversity}
\description{Compute true diversities of individual samples from a vector of relative abundances or a matrix (OTU table) containing multiple samples. If a tree is provided, the computed true diversities consider the phylogenetic relations between OTUs.}
\usage{true.div(abund,qvalue,tree)}
\arguments{
  \item{abund}{A vector or a matrix indicating the relative abundances of one or multiple samples. If a matrix is provided, columns should be samples and rows OTUs. The relative abundances of each sample (vector or each column of the matrix) must sum to 1.}
  \item{qvalue}{A positive number (>0). It can be an integer or contain decimals.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
}
\details{True diversities, Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs.}
\examples{
true.div(abund=vector,qvalue=1)
true.div(abund=vector,qvalue=1,tree=tree)
true.div(vector,1,tree)
true.div(otu.table,0,tree)
}
