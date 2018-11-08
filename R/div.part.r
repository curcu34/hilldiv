div.part <- function(otutable,qvalue,hierarchy) {
  
#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}
if(missing(hierarchy)) warning("Assuming a two-level hierarchy: 1) sample, 2) total dataset")

#Function for 2-level hierarchy
if(missing(hierarchy)){
alpha <- DiverHill::alpha.div(otutable,qvalue)
gamma <- DiverHill::gamma.div(otutable,qvalue)
beta <- gamma/alpha
N <- ncol(otutable)
homogeneity <- ((1/beta) - 1/N)/(1-1/N)
overlap <-((1/beta)^(1-qvalue) - (1/N)^(1-qvalue)) / (1 - (1/N)^(1-qvalue))
turnover <- (beta - 1)/(N-1)
results <- list("Sample_size" = N, "Alpha_diversity" = alpha, "Gamma_diversity" = gamma, "Beta_diversity" = beta, "Homogeneity" = homogeneity, "Overlap" = overlap, "Turnover" = turnover)
return(results)
}
  
#Function for 3-level hierarchy
if(ncol(hierarchy) == 2){

results <- c("Allright")
return(results)
}
