\name{beta.dis}
\alias{beta.dis}
\title{(Dis)similarity measures derived from Hill number beta diversities}
\usage{beta.dis(beta,qvalue,N,metric,type)}
\arguments{
  \item{beta}{A numeric beta value or a 'div.part' function output object.}
  \item{qvalue}{A positive number (>0). It can be an integer or contain decimals. The argument is ovewritten if a 'div.part' object is used.}
  \item{N}{Integer. Sample size, the number of sampling units to be used to compute the (dis)similarity measures. The argument is ovewritten if a 'div.part' object is used}.
  \item{metric}{A character object containing "C", "U", "V" or "S". See the explanation of each metric below. metric="U".}
  \item{type}{A character object containing either "similarity" or "dissimilarity". Defines the 'orientation' of the metric. type="dissimilarity"}
}
\details{It is possible to compute multiple (dis)similarity measurements derived from beta diversities, both for Hill numbers and phylogenetic Hill numbers. Four classes of similarity measures derived from Hill numbers beta diversities have been developed. The Jaccard-type overlap (UqN) quantifies the effective proportion of OTUs or lineages in a system that are shared across all subsystems. The Sørensen-type turnover-complement (VqN) is the complement of the Sørensen-type turnover, which quantifies the normalized OTU turnover rate with respect to the average subsystem (i.e. alpha). Finally, the Jaccard-type turnover-complement (SqN) is the complement of the Jaccard-type turnover, which quantifies the normalized OTU turnover rate with respect to the whole system (i.e. gamma). Dissimilarity measures can be obtained by calculating their one-complements (1-XqN).}
\examples{
beta.dis(beta=4.5,qvalue=1,N=8)
beta.dis(beta=4.5,qvalue=1,N=8,metric="C",type="similarity")
div.part.object  <- div.part(otu.table,qvalue=0,tree=tree)
beta.dis(div.part.object)
beta.dis(div.part.object,metric="S",type="similarity")
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
Chao, A., Chiu, C.‐H., & Hsieh, T. C. (2012). Proposing a resolution to de‐ bates on diversity partitioning. Ecology, 93, 2037–2051\cr\cr
Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439.
}
