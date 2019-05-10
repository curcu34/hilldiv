
\name{hill.div}
\alias{hill.div}
\title{True diversity}
\description{Compute Hill numbers of individual samples from a vector of relative abundances or a matrix (OTU table) containing multiple samples. If an ultrametic OTU tree is provided, phylogenetic Hill numbers are computed.}
\usage{true.div(abund,qvalue,tree,type)}
\arguments{
  \item{abund}{A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns refer to samples and rows to OTUs.}
  \item{qvalue}{A positive number (>=0). It can be an integer or contain decimals.}
  \item{tree}{An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{type}{Either "abundance" for abundance-based Hill number computation or "incidence" for incidence-based Hill number computation. In the latter case, the abund object needs to be an OTU table, as incidence-based Hill number computation requires multiple samples. Default: type="abundance".}
}
\details{Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs. If the relative abundances of each sample (vector or each column of the matrix) do not sum to 1, TSS-normalisation is conducted.}
\examples{
hill.div(otu.vector,0)
hill.div(otu.table,1)
hill.div(otu.table,qvalue=2)
hill.div(otu.table,1,tree)
hill.div(otu.table,1,tree,type="incidence")
hill.div(otu.table,qvalue=2,tree=tree)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
Jost, L. (2006). Entropy and diversity. Oikos, 113, 363–375.\cr\cr
Hill, M. O. (1973). Diversity and evenness: a unifying notation and its con‐ sequences. Ecology, 54, 427–432.
}
