
\name{div.test.plot}
\alias{div.test.plot}
\title{Plot of diversity comparison between groups of samples}
\usage{div.test.plot(divtest,chart,colour)}
\arguments{
  \item{divtest}{Object outputed by the div.test() function}
  \item{chart}{Chart type, either 'box' for boxplot, 'jitter' for jitter plot or 'violin' for violin plot. chart="box"}
  \item{colour}{Vector of RGB colours | The number of vector items (colours, e.g. '#34k235'), must equal the number of groups that are intended to plot.}
}
  \details{Visualises the diversity differences between groups of samples based on the information contained by the object outputed by the div.test() function.}
\examples{
contrast.div.q0 <- div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
div.test.plot(contrast.div.q0,chart="jitter")
div.test.plot(contrast.div.q0,chart="violin")
}
