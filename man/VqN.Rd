\name{VqN}
\alias{VqN}
\title{VqN - Sørensen-type turnover-complement}
\usage{VqN(beta,N)}
\arguments{
  \item{beta}{A beta diversity value based on Hill numbers.}
  \item{N}{Integer | Sample size.}
}
\details{The Sørensen-type turnover-complement is the complement of the Sørensen‐type turnover, which quantifies the normalized OTU turnover rate with respect to the average subsystem (i.e., alpha), thus provides the proportion of a typical subsystem that changes across subsystems. VqN is integrated in the functions beta.dis() and pair.dis().}
\examples{
VqN(beta=1.24,N=2)
VqN(1.24,2)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.
Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439. 
Harrison, S., Ross, S. J., & Lawton, J. H. (1992). Beta diversity on geo‐ graphic gradients in Britain. Journal of Animal Ecology, 61, 151–158.
}
