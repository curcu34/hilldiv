\name{CqN}
\alias{CqN}
\title{CqN - Sørensen-type overlap}
\usage{CqN(beta,qvalue,N)}
\usage{div.part(otutable,qvalue,hierarchy,tree,type)}
\arguments{
  \item{beta}{A beta diversity value based on Hill numbers.}
  \item{qvalue}{A positive number (>=0). It can be an integer or contain decimals.}
  \item{N}{Integer | Sample size.}
}
\details{The Sørensen-type overlap quantifies the effective average proportion of a sub‐systems OTUs (or lineages in the case of phylodiversities) that is shared across all subsystems. This is thus a metric that quantifies overlap from the subsystems perspective. Its corresponding dissimilarity measure (1 − CqN) quantifies the effective average proportion of nonshared OTUs or lineages in a system. CqN is integrated in the functions beta.dis() and pair.dis().}
\examples{
CqN(beta=1.24,qvalue=1,N=2)
CqN(1.24,1,2)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
Chao, A., Chiu, C.‐H., & Hsieh, T.C. (2012). Proposing a resolution to de‐ bates on diversity partitioning. Ecology, 93, 2037–2051.\cr\cr
Chao, A., Jost, L., Chiang, S.C., Jiang, Y.‐H., & Chazdon, R.L. (2008). A two‐stage probabilistic approach to multiple‐commu‐ nity similarity indices. Biometrics, 64, 1178–1186.
}
