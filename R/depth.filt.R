#' Sequencing depth filtering
#' @title Sequencing depth filtering
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords coverage sequencing
#' @description Filter samples based on a minimum sequencing depth.
#' @param otutable An OTU table (matrix/data.frame) indicating the absolute OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs.
#' @param threshold A number indicating the minimum sequencing depth required to keep the sample.
#' @seealso \code{\link{depth.cov}}, \code{\link{copy.filt}}
#' @examples
#' depth.filt(otu.table,5000)
#' depth.filt(otu.table,threshold=20000)
#' @references
#' Alberdi A, Aizpurua O, Bohmann K, Gopalakrishnan S, Lynggaard C, Nielsen M, Gilbert MTP. 2019. Promises and pitfalls of using high-throughput sequencing for diet analysis. Molecular Ecology Resources, 19(2), 327-348.\cr\cr
#' @export

depth.filt <- function(otutable,threshold){
if(colSums(otutable)[1] == 1)stop("OTU table contains relative abundances. This function is only meaningful for absolute abundances.")
otutable.filt <- otutable[,colSums(otutable)>threshold]
return(otutable.filt)
}
