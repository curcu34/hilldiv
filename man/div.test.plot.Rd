
\name{div.test.plot}
\alias{div.test.plot}
\title{Plot of diversity comparison between groups of samples}
\usage{div.test.plot(divtest,chart,colour)}
\arguments{
  \item{divtest}{Object outputed by the div.test() function}
  \item{chart}{Chart type, either 'box' for boxplot, 'jitter' for jitter plot or 'violin' for violin plot. chart="box"}
  \item{colour}{Vector of RGB colours | The number of vector items (colours, e.g. '#34k235'), must equal the number of groups that are intended to plot.}
  \item{stat}{Logical (TRUE/FALSE) | If 'TRUE' pairwise mean comparison significance values will be ploted.}
  \item{comb}{List of pairwise combinations. If missing all combinations will be plotted. e.g. list(c("Myotis myotis","Myotis capaccinii"),c("Myotis myotis","Myotis daubentonii"))}
  \item{symbol}{Logical (TRUE/FALSE) | If 'TRUE' symbols rather than p-values will be shown. ns: p > 0.05; *: p <= 0.05; **: p <= 0.01; ***: p <= 0.001; ****: p <= 0.0001.}
  \item{flip}{Logical (TRUE/FALSE) | If 'TRUE' the chart will be flipped 90 degrees.}
}
  \details{Visualises the diversity differences between groups of samples based on the information contained by the object outputed by the div.test() function.}
\examples{
contrast.div.q0 <- div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
div.test.plot(contrast.div.q0,chart="jitter")
div.test.plot(contrast.div.q0,chart="violin")
div.test.plot(div.test.result,chart="jitter",stat=TRUE,flip=TRUE)
div.test.plot(div.test.result,stat=TRUE,comb=list(c("Myotis myotis","Myotis capaccinii")),symbol=TRUE)
}
