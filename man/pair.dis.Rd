\name{pair.dis}
\alias{pair.dis}
\title{Pairwise dissimilarity}
\usage{pair.dis(otutable,qvalue,tree,hierarchy,level,measure)}
\arguments{
  \item{otutable}{A matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns refer to samples and rows to OTUs.}
  \item{qvalue}{A positive number (>=0). It can be an integer or contain decimals.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{hierarchy}{2-column matrix | The first column lists the sample names while the second lists the groups. If provided, group profiles are plotted instead of individual profiles.}
  \item{level}{Either 1 or 2 | If '1' dissimilarities are computed across samples, while if '2' dissimilarities are computed across groups, as specified by the hierarchy table. level=1}
}
\examples{
pair.dis(otu.table,qvalue=1)
pair.dis(otu.table,qvalue=1,tree=tree)
pair.dis(otu.table,qvalue=0,hierarchy=hierarchy.table,level="2")
}

