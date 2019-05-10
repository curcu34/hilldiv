\name{SqN}
\alias{SqN}
\title{SqN - Jaccard-type turnover-complement }
\usage{SqN(beta,N)}
\arguments{
  \item{beta}{A beta diversity value based on Hill numbers.}
  \item{N}{Integer | Sample size.}
}
\details{The Jaccard-type turnover-complement is thecomplement of the Jaccard‐type turnover, which quantifies the normalized OTU turnover rate with respect to the whole system (i.e. gamma). SqN is integrated in the functions beta.dis() and pair.dis().}
\examples{
SqN(beta=1.24,N=2)
SqN(1.24,2)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439.\cr\cr
Harrison, S., Ross, S. J., & Lawton, J. H. (1992). Beta diversity on geo‐ graphic gradients in Britain. Journal of Animal Ecology, 61, 151–158.
}
