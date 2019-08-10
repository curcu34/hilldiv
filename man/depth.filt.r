\name{depth.filt}
\alias{depth.filt}
\title{Filter samples based on a minimum sequencing depth}
\usage{depth.filt(otutable,threshold)}
\arguments{
  \item{otutable}{An OTU table (matrix/data.frame) indicating the absolute OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs.}
  \item{threshold}{A number indicating the minimum sequencing depth required to keep the sample.}
}
\details{Metabarcoding studies usually require a minimum sequencing depth is usually required as a quality control.}
\examples{
depth.filt(otu.table,5000)
depth.filt(otu.table,threshold=20000)
}
\references{
Alberdi A, Aizpurua O, Bohmann K, Gopalakrishnan S, Lynggaard C, Nielsen M, Gilbert MTP. 2019. Promises and pitfalls of using high-throughput sequencing for diet analysis. Molecular Ecology Resources, 19(2), 327-348.\cr\cr
}
