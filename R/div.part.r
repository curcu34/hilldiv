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
if(identical(sort(rownames(otutable)),sort(hierarchy[,1])) == FALSE) stop("OTU names in the OTU table and the hierarchy table do not match")
colnames(hierarchy) <- c("L1","L2")
  
L1_div <- DiverHill::alpha.div(otutable,qvalue)
L3_div <- DiverHill::gamma.div(otutable,qvalue)

otutable.L2 <- merge(t(otutable),hierarchy, by.x="row.names",by.y="L1")
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- aggregate(otutable.L2[,-ncol(otutable.L2)], by=list(otutable.L2[,ncol(otutable.L2)]), FUN=sum)
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- t(sweep(otutable.L2, 1, rowSums(otutable.L2), FUN="/"))
L2_div <- DiverHill::alpha.div(otutable.L2,qvalue)

beta_1 <- L2_div/L1_div
N1 <- ncol(otutable)
homogeneity_1 <- ((1/beta_1) - 1/N1)/(1-1/N1)
overlap_1 <-((1/beta_1)^(1-qvalue) - (1/N1)^(1-qvalue)) / (1 - (1/N1)^(1-qvalue))
turnover_1 <- (beta_1 - 1)/(N1-1)
  
beta_2 <- L3_div/L2_div
N2 <- ncol(otutable.L2)
homogeneity_2 <- ((1/beta_2) - 1/N2)/(1-1/N2)
overlap_2 <-((1/beta_2)^(1-qvalue) - (1/N2)^(1-qvalue)) / (1 - (1/N2)^(1-qvalue))
turnover_2 <- (beta_2 - 1)/(N2-1)
  
results <- list("L1_size" = N1, "L2_size" = N2, "L1_diversity" = L1_div, "L2_diversity" = L2_div, "L3_diversity" = L3_div, "Beta_diversity_L1-2" = beta_1, "Beta_diversity_L2-3" = beta_2,  "Homogeneity_L1-2" = homogeneity_1, "Homogeneity_L2-3" = homogeneity_2,"Overlap_L1-2" = overlap_1, "Overlap_L2-3" = overlap_2,"Turnover_L1-2" = turnover_1, "Turnover_L2-3" = turnover_2)
return(results)
}
  
#Function for 4-level hierarchy
if(ncol(hierarchy) == 3){

results <- list("L1_size" = N1, "L2_size" = N2, "L1_diversity" = L1_div, "L2_diversity" = L2_div, "L3_diversity" = L3_div, "Beta_diversity_L1-2" = beta_1, "Beta_diversity_L2-3" = beta_2,  "Homogeneity_L1-2" = homogeneity_1, "Homogeneity_L2-3" = homogeneity_2,"Overlap_L1-2" = overlap_1, "Overlap_L2-3" = overlap_2,"Turnover_L1-2" = turnover_1, "Turnover_L2-3" = turnover_2)
return(results)
} 
  
#Error if more hierarchical levels
if(ncol(hierarchy) > 3){
stop("The maximum number of hierarchical levels allowed is 4")
} 

}
