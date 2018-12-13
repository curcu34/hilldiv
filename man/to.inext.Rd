\name{to.inext}
\alias{to.inext}
\title{Transform OTU table to iNext 'incidence_raw' format}
\usage{}
\arguments{}
\examples{
sp.lists <- to.inext(otu.table.counts,sample.species)
sp.inext <- iNEXT(sp.lists, q=0, datatype="incidence_freq")
ggiNEXT(sp.inext, type=1, facet.var="order", color.var="sp")
}
