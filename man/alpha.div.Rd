\name{alpha.div}
\alias{alpha.div}
\title{Alpha diversity}
\description{Compute alpha diversity of a system from a matrix (OTU table) containing multiple samples. If a tree is provided, the computed alpha diversity consider the phylogenetic relations between OTUs.}
\usage{alpha.div(otutable,qvalue,weight)}
\arguments{
  \item{otutable}{A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs. The relative abundances of each sample (columns) must sum to 1.}
  \item{qvalue}{A positive number (>0). It can be an integer or contain decimals.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
}
\details{Alpha diversity of a system from a matrix (OTU table) containing multiple samples.}
\examples{
alpha.div(otutable=otu.table,qvalue=1,weight=weight.vector)
alpha.div(otutable=otu.table,qvalue=1,tree=tree)
alpha.div(otu.table,1,eight.vector,tree)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
Chao, A., Chiu, C.‐H., & Hsieh, T. C. (2012). Proposing a resolution to de‐ bates on diversity partitioning. Ecology, 93, 2037–2051\cr\cr
Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439.
}
