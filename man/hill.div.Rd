
\name{hill.div}
\alias{hill.div}
\title{True diversity}
\description{Compute true diversities of individual samples from a vector of relative abundances or a matrix (OTU table) containing multiple samples. If a tree is provided, the computed true diversities consider the phylogenetic relations between OTUs.}
\usage{true.div(abund,qvalue,tree)}
\arguments{
  \item{abund}{A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns refer to samples and rows to OTUs.}
  \item{qvalue}{A positive number (>=0). It can be an integer or contain decimals.}
  \item{tree}{An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{type}{Either "abundance" for abundance-based Hill number computation or "incidence" for incidence-based Hill number computation. In the latter case, the abund object needs to be a matrix/data.frame, as incidence-based Hill number computation requires multiple samples. Default: type="abundance".}
}
\details{True diversities, Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs. If the relative abundances of each sample (vector or each column of the matrix) do not sum to 1, TSS-normalisation is conducted.}
\examples{
true.div(abund=vector,qvalue=1)
true.div(abund=vector,qvalue=1,tree=tree)
true.div(vector,1,tree)
true.div(otu.table,0,tree)
}
