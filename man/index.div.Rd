
\name{index.div}
\alias{index.div}
\title{Diversity indices}
\description{Computes common diversity indices related to Hill numbers. If the input is a vector, the function computes the indices of a single sample, while if the input is a matrix (OTU table), the function computes individual diversity indices for each sample (column). An ultrametic OTU tree is required for computing phylogenetic diversity indices (Faith's PD, Allen's H and Rao's Q). If the relative abundances of each sample (vector or each column of the matrix) do not sum to 1, TSS-normalisation is appplied.}
\usage{hill.div(abund,tree,index)}
\arguments{
  \item{abund}{A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.}
  \item{tree}{An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{index}{Diversity index to be computed ("richness", "shannon", "simpson", "faith", "allen", "rao"). Default without tree argument: index="richness". Default with tree argument: index="faith".}
}
\details{The six diversity indices the index.div() function computes are closely related to Hill numbers.}
\examples{
index.div(otu.vector)
index.div(otu.table,index="shannon")
index.div(otu.table,tree=bat.diet.tree)
index.div(bat.diet.otutable,tree=bat.diet.tree,index="allen")
index.div(bat.diet.otutable,bat.diet.tree,"rao")
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
Jost, L. (2006). Entropy and diversity. Oikos, 113, 363–375.\cr\cr
Rao, C. R. (1982). Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology, 21, 24–43.\cr\cr
Shannon, C. E. (1948). A mathematical theory of communication. The Bell System Technical Journal, 27, 379–423.\cr\cr
}
