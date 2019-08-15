#' Diversity partitioning (based on Hill numbers)
#' @title Diversity partitioning
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords alpha gamma beta hill
#' @description Diversity partitioning based on Hill numbers. the Hill numbers framework enables the diversity of a system to be partitioned following the multiplicative definition (Jost, 2007; Chao, Chiu, & Hsieh, 2012). Alpha diversity is obtained by computing the Hill numbers from the averaged basic sums of the samples, while gamma diversity is obtained by taking the average of OTU relative abundances across samples, and then computing the Hill numbers of the pooled system. The division of gamma diversity by alpha diversity yields the beta diversity, which quantifies how many times richer an entire system is in effective OTUs (gamma diversity) than its constituent samples are on average (alpha diversity). However, the Hill numbers beta diversity can also be considered an actual diversity value, as the same metric also measures the effective number of equally large and completely distinct samples in a system.
#' @param otutable An OTU table (matrix/data.frame) indicating the absolute or relative OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs.
#' @param qvalue A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals.
#' @param hierarchy A two-column matrix indicating the relation between samples (first column) and groups (second column).
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match.
#' @param type Either "abundance" for abundance-based Hill number computation or "incidence" for incidence-based Hill number computation. Incidence-based Hill number computation requires a hierarchy table.
#' @return A list object containing general information about the diversity partitioning and L1 (level 1 = alpha), L2 (level 2 = gamma) and beta diversities.
#' @seealso \code{\link{div.part}}, \code{\link{gamma.div}}, \code{\link{match.data}}
#' @examples
#' div.part(otu.table,qvalue=1)
#' div.part(otu.table,qvalue=1,type="abundance")
#' div.part(otu.table,qvalue=0,tree=tree)
#' div.part(otu.table,qvalue=0,type="incidence",hierarchy=hierarchy.table)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' Chao, A., Chiu, C.‐H., & Hsieh, T. C. (2012). Proposing a resolution to de‐ bates on diversity partitioning. Ecology, 93, 2037–2051
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439.
#' @export

div.part <- function(otutable,qvalue,tree,type,hierarchy) {

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
###########
#Function for 2-level hierarchy
###########

#L1 and L2 diversities
if(missing(hierarchy)){
if(type == "incidence") stop("Diversity partitioning based on incidence data requires a hierarchy table")
if(missing(tree)){
  L1_div <- hilldiv::alpha.div(otutable,qvalue)
  L2_div <- hilldiv::gamma.div(otutable,qvalue)
  }else{
  L1_div <- hilldiv::alpha.div(otutable,qvalue,tree)
  L2_div <- hilldiv::gamma.div(otutable,qvalue,tree)
}

#Beta and sample size
beta <- L2_div/L1_div
N <- ncol(otutable)

#Return values
if (qvalue==0.99999) {qvalue=1}
results <- list("Hierarchical_levels" = 2,"Type" = type,"Order_diversity" = qvalue,"Sample_size" = N, "L1_diversity" = L1_div, "L2_diversity" = L2_div, "Beta_diversity" = beta)
return(results)
}

#Pre-processing
if(ncol(hierarchy) == 2){
hierarchy[,1] <- as.character(hierarchy[,1])
hierarchy[,2] <- as.character(hierarchy[,2])
if(identical(sort(colnames(otutable)),sort(hierarchy[,1])) == FALSE) stop("OTU names in the OTU table and the hierarchy table do not match")
colnames(hierarchy) <- c("L1","L2")

if(type == "abundance"){
#####
# Abundance-based
#####

  #L2 data processing
  otutable.L2 <- merge(t(otutable),hierarchy, by.x="row.names",by.y="L1")
    rownames(otutable.L2) <- otutable.L2[,1]
  otutable.L2 <- otutable.L2[,-1]
  otutable.L2 <- aggregate(otutable.L2[,-ncol(otutable.L2)], by=list(otutable.L2[,ncol(otutable.L2)]), FUN=sum)
  rownames(otutable.L2) <- otutable.L2[,1]
  otutable.L2 <- otutable.L2[,-1]
  otutable.L2 <- t(sweep(otutable.L2, 1, rowSums(otutable.L2), FUN="/"))

  #Get L2 weights
  weights.L2 <- table(hierarchy[,2])/sum(table(hierarchy[,2]))

  #L1 and L3 diversities
  if(missing(tree)){
    L1_div <- hilldiv::alpha.div(otutable,qvalue)
    L3_div <- hilldiv::gamma.div(otutable,qvalue)
    }else{
   L1_div <- hilldiv::alpha.div(otutable,qvalue,tree)
   L3_div <- hilldiv::gamma.div(otutable,qvalue,tree)
  }

  #L2 diversity
  if(missing(tree)){
    L2_div <- hilldiv::alpha.div(otutable.L2,qvalue,weight=weights.L2)
    }else{
    L2_div <- hilldiv::alpha.div(otutable.L2,qvalue,tree,weight=weights.L2)
  }

  #L1-L2, beta and sample size
  beta_1 <- L2_div/L1_div
  N1 <- ncol(otutable)

  #L2-L3, beta and sample size
  beta_2 <- L3_div/L2_div
  N2 <- ncol(otutable.L2)

  #Return values
  if (qvalue==0.99999) {qvalue=1}
  results <- list("Hierarchical_levels" = 3,"Type" = type,"Order_diversity" = qvalue, "L1_Sample_size" = N1, "L2_Sample_size" = N2, "L1_diversity" = L1_div, "L2_diversity" = L2_div, "L3_diversity" = L3_div, "Beta_diversity_L1_2" = beta_1, "Beta_diversity_L2_3" = beta_2)
  return(results)

}else if(type == "incidence"){
#####
# Incidence-based
#####

  #L2 data processing
  otutable.L2 <- merge(t(otutable),hierarchy, by.x="row.names",by.y="L1")
  rownames(otutable.L2) <- otutable.L2[,1]
  otutable.L2 <- otutable.L2[,-1]
  otutable.L2 <- aggregate(otutable.L2[,-ncol(otutable.L2)], by=list(otutable.L2[,ncol(otutable.L2)]), function(x) sum(x != 0))
  rownames(otutable.L2) <- otutable.L2[,1]
  otutable.L2 <- otutable.L2[,-1]
  otutable.L2 <- t(sweep(otutable.L2, 1, rowSums(otutable.L2), FUN="/"))

  #Get L2 weights
  weights.L2 <- table(hierarchy[,2])/sum(table(hierarchy[,2]))

  #L2 and L3 diversities
  if(missing(tree)){
    L2_div <- hilldiv::alpha.div(otutable.L2,qvalue,weight=weights.L2)
    L3_div <- hilldiv::gamma.div(otutable.L2,qvalue)
    }else{
    L2_div <- hilldiv::alpha.div(otutable.L2,qvalue,tree,weight=weights.L2)
    L3_div <- hilldiv::gamma.div(otutable.L2,qvalue,tree)
  }

  #Beta and sample size
  beta <- L3_div/L2_div
  N <- length(weights.L2)

  #Return values
  if (qvalue==0.99999) {qvalue=1}
  results <- list("Hierarchical_levels" = 2,"Type" = type,"Order_diversity" = qvalue,"Sample_size" = N, "L2_diversity" = L2_div, "L3_diversity" = L3_div, "Beta_diversity" = beta)
  return(results)
}
}
}
