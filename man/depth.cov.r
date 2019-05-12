\name{depth.cov}
\alias{depth.cov}
\title{Coverage of the estimated Hill numbers at different orders of diversity}
\usage{depth.cov(abund,qvalue)}
\arguments{
  \item{abund}{A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.}
  \item{qvalue}{A positive number (>=0). It can be an integer or contain decimals.}
}
\details{Assessment of the sequencing depth per sample. Computes both observed and estimated (based on Chao & Jost 2015) Hill numbers, and calculates the percentage of estimated diversity coverage per sample.}
\examples{
depth.cov(otu.table,0)
depth.cov(otu.table,qvalue=1)
}
\references{
Chao, A. & Jost, L. (2015) Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution, 6, 873-882.\cr\cr
Jost, L. (2006). Entropy and diversity. Oikos, 113, 363–375.\cr\cr
Hill, M. O. (1973). Diversity and evenness: a unifying notation and its con‐ sequences. Ecology, 54, 427–432.
}
