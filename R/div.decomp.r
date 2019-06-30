div.decomp <- function(otutable,qvalue,tree,type,hierarchy) {

#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(qvalue==1) {qvalue=0.99999}
if(missing(type)){type = "abundance"}
if(type == "incidence"){
  if(missing(hierarchy)) stop("Diversity partitioning based on incidence data requires a hierarchy table")
}


}
