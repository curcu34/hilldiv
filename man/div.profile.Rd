\name{div.profile}
\alias{div.profile}
\title{Diversity profile}
\usage{div.profile(abund,tree,order,value)}
\arguments{
  \item{abund}{A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs. The relative abundances of each sample (columns) must sum to 1.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{order}{A vector of sequential orders of diversity (usually from 0 to 5). order=seq(from = 0, to = 5, by = (0.1))}
  \item{values}{Logical (TRUE/FALSE) | If 'TRUE' the Hill numbers for each order of diversity for each sample or group (if hierarchy table is provided) are printed on screen. values="FALSE"}
  \item{hierarchy}{2-column matrix | The first column lists the sample names while the second lists the groups. If provided, group profiles are plotted instead of individual profiles.}
  \item{level}{Either "alpha" or "gamma" | Only if hierachy is used. If 'alpha', alpha diversity of the groups is plotted instead of gamma diversity. level="gamma"}
  \item{log}{Logical (TRUE/FALSE) | If 'TRUE' the Hill numbers are transformed to the logarithmic scale. This is useful when there are large differences between q values (e.g. sharp drop from q=0 to q=1), which might complicate visualization. log="FALSE"}
}
\details{The diversity profile enables visualizing the diversity distribution of a sample or a group of samples at a glance.}
\examples{
div.profile(abund=otu.table[,1],tree=tree,order=seq(from = 0, to = 5, by = (0.1)),value="TRUE")
div.profile(vector)
div.profile(otu.table,tree)
}
