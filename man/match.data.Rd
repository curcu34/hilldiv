
\name{match.data}
\alias{match.data}
\title{Match data}
\description{Filter OTU tables and OTU trees to match OTUs present in both data files.}
\usage{match.data(otutable,tree)}
\arguments{
  \item{otutable}{An OTU table (matrix/data.frame) indicating the absolute OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs.}
  \item{tree}{An ultrametic tree of class 'phylo'.}
  \item{output}{Whether to output a filtered OTU table (matrix) or a filtered OTU tree (phylo).}
}
\details{Auxiliary function to match OTUs present in OTU tables and OTU trees and ensure phylogenetic Hill numbers are correctly computed.}
\examples{
match.data(bat.diet.otutable,bat.diet.tree,output="otutable")
match.data(bat.diet.otutable,bat.diet.tree,output="tree")
}
