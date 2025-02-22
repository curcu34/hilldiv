% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depth.cov.R
\name{depth.cov}
\alias{depth.cov}
\title{Depth coverage assessment}
\usage{
depth.cov(abund, qvalue)
}
\arguments{
\item{abund}{A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.}

\item{qvalue}{A positive integer or decimal number (>=0), usually between 0 and 3.}
}
\value{
A matrix with observed diversity,  estimated diversities and coverage
}
\description{
Coverage of the estimated Hill numbers at different orders of diversity.
}
\details{
Depth coverage assessment
}
\examples{
depth.cov(otu.table,0)
depth.cov(otu.table,qvalue=1)
}
\references{
Chao, A. & Jost, L. (2015) Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution, 6, 873-882.\cr\cr
Jost, L. (2006). Entropy and diversity. Oikos, 113, 363–375.
Hill, M. O. (1973). Diversity and evenness: a unifying notation and its con‐ sequences. Ecology, 54, 427–432.
}
\seealso{
\code{\link{hill.div}}, \code{\link{depth.filt}}
}
\author{
Antton Alberdi, \email{anttonalberdi@gmail.com}
}
\keyword{coverage}
\keyword{sequencing}
