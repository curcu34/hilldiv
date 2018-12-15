\name{hill.intext}
\alias{hill.intext}
\title{Interpolation and extrapolation of Hill numbers and phylogenetic Hill numbers}
\description{The function wraps the R packages iNEXT and iNextPD to carry out interpolation and extrapolation of Hill numbers and phylogenetic Hill numbers.}
\usage{alpha.div(otutable,qvalue,weight)}
\arguments{
  \item{otutable}{A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs. The relative abundances of each sample (columns) must sum to 1.}
  \item{qvalue}{Integers 0, 1 or 2 | qvalue=1}
  \item{tree}{An ultrametric tree of class 'phylo' or 'phylog' | The tip labels must match the row names of the OTU table. In large phylogenetic trees the transformation from 'phylo' to 'phylog' can take a long time, so it is advised to transform the tree using the function 'phylo.to.phylog' beforehand.}
  \item{hierarchy}{2-column matrix | The first column lists the sample names while the second lists the groups. If provided, group profiles are plotted instead of individual profiles.}
  \item{type}{Either "incidence" or "abundance" | If 'incidence', system-level or group-level (if hierarchy file is provided) incidence-based rarefaction and extrapolation is conducted, while if 'abundance', subsystem/sample-level rarefaction and extrapolation is conducted. type="incidence"}
  \item{output}{'diversity', 'completeness', 'report' | String indicating the type of results to output. 'diversity' yields a diversity accumulation curve plot, 'completeness' yields a completeness curve plot, 'report' outputs the iNEXT/iNextPD object. output="diversity"}
  \item{colour}{Vector of RGB codes | Vector containing specific colour codes to be used. colour=c("#324353","#745431","#266576")}
  \item{size}{Vector of integers | Vector containing increasing sample sizes. If not used, the Hill numbers of all samples/groups are extrapolated three times te size of the largest sample in 20 steps. size=seq(1,maxsize*3,round(maxsize*3/20))}
}
\examples{
hill.intext(otutable=otu.table,qvalue=1,tree=otu.tree,hierarchy=hierarchy.table,output="completeness",size=seq(1,maxsize*3,round(maxsize*3/20)))
}
