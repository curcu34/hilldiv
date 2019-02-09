\name{pair.dis.plot}
\alias{pair.dis.plot}
\title{Visualisation of pairwise dissimilarities}
\usage{pair.dis.plot(distance,hierarchy,type,level,colour,magnify)}
\arguments{
  \item{distance}{Matrix of pairwise dissimilarities, usually one of the matrices listed in the output object of the pair.dis() function}
  \item{hierarchy}{2-column matrix | The first column lists the sample names while the second lists the groups. If provided, group profiles are plotted instead of individual profiles.}
  \item{type}{Whether to plot a NMDS or qgraph chart. type="NMDS"}
  \item{level}{Either 1 or 2 | If '1', samples are not grouped (each sample is a different colour), while if '2', samples are coloured by groups. level=1}
  \item{colour}{Vector of RGB colours | The number of vector items (colours, e.g. '#34k235'), must equal the number of samples or groups that are intended to plot.}
  \item{magnify}{Only relevant for qgraph. Whether the pairwise dissimilarity values are transformed to 0-1 scale, 0 corresponding to the minimum dissimilarity and 1 to the maximum dissimilarity value. magnify=FALSE}
}
\examples{
pair.div.q0.L2 <- pair.dis(otu.table,qvalue=0,hierarchy=hierarchy.table,level="2")
pair.dis.plot(pair.div.q0.L2$L2_CqN,hierarchy=hierarchy.table,type="NMDS",level=2)
pair.dis.plot(pair.div.q0.L2$L2_CqN,hierarchy=hierarchy.table,type="qgraph",level=2,magnify=TRUE)
}
