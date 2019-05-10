
\name{div.part}
\alias{div.part}
\title{Diversity partitioning}
\usage{div.part(otutable,qvalue,hierarchy,tree,type)}
\arguments{
  \item{otutable}{A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs.}
  \item{qvalue}{A positive number (>=0). It can be an integer or contain decimals.}
  \item{hierarchy}{2-column matrix | The first column lists the sample names while the second lists the groups. If provided, group profiles are plotted instead of individual profiles.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{type}{Either "abundance" for abundance-based Hill number computation or "incidence" for incidence-based Hill number computation. Incidence-based Hill number computation requires a hierarchy table. Default: type="abundance".}
}
\details{the Hill numbers framework enables the diversity of a system to be partitioned following the multiplicative definition (Jost, 2007; Chao, Chiu, & Hsieh, 2012). Alpha diversity is obtained by computing the Hill numbers from the averaged basic sums of the samples, while gamma diversity is obtained by taking the average of OTU relative abundances across samples, and then computing the Hill numbers of the pooled system. The division of gamma diversity by alpha diversity yields the beta diversity, which quantifies how many times richer an entire system is in effective OTUs (gamma diversity) than its constituent samples are on average (alpha diversity). However, the Hill numbers beta diversity can also be considered an actual diversity value, as the same metric also measures the effective number of equally large and completely distinct samples in a system.}
\examples{
div.part(otu.table,qvalue=1)
div.part(otu.table,qvalue=1,type="abundance")
div.part(otu.table,qvalue=0,tree=tree)
div.part(otu.table,qvalue=0,type="incidence",hierarchy=hierarchy.table)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
Chao, A., Chiu, C.‐H., & Hsieh, T. C. (2012). Proposing a resolution to de‐ bates on diversity partitioning. Ecology, 93, 2037–2051\cr\cr
Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439.
}
