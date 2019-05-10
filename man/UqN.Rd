\name{UqN}
\alias{UqN}
\title{UqN - Jaccard-type overlap}
\usage{UqN(beta,qvalue,N)}
\arguments{
  \item{beta}{A beta diversity value based on Hill numbers.}
  \item{qvalue}{A positive number (>=0). It can be an integer or contain decimals.}
  \item{N}{Integer | Sample size.}
}
\details{The Jaccard-type overlap quantifies the effective proportion of OTUs or lineages in a system that are shared across all subsystems. Hence, this metric quantifies overlap from the perspective of the overall system. Its corresponding dissimilarity (1 − UqN) quantifies the effective proportion of nonshared OTUs or lineages in the overall system. UqN is integrated in the functions beta.dis() and pair.dis().}
\examples{
UqN(beta=1.24,qvalue=1,N=2)
UqN(1.24,1,2)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.
Chao, A., Chiu, C.‐H., & Hsieh, T.C. (2012). Proposing a resolution to de‐ bates on diversity partitioning. Ecology, 93, 2037–2051.
Chao, A., Jost, L., Chiang, S.C., Jiang, Y.‐H., & Chazdon, R.L. (2008). A two‐stage probabilistic approach to multiple‐commu‐ nity similarity indices. Biometrics, 64, 1178–1186.
}
