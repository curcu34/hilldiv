
\name{tree.depth}
\alias{tree.depth}
\title{Tree depth}
\description{Computes phylogenetic tree depth based from a phylogenetic tree and a vector of (relative) abundances.}
\usage{hill.div(tree,vector)}
\arguments{
  \item{tree}{An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{vector}{A vector indicating the absolute or relative OTU abundances of one sample.}
}
\details{Auxiliary function to compute tree depths.}
\examples{
tree.depth(bat.diet.tree,otu.vector)
tree.depth(tree=bat.diet.tree,vector=otu.vector)
}
\references{
Jost, L. (2006). Entropy and diversity. Oikos, 113, 363–375.\cr\cr
}
