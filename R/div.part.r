div.part <- function(otutable,qvalue,hierarchy,tree) {
  
#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}

###########
#Function for 2-level hierarchy
########### 

#L1 and L2 diversities
if(missing(hierarchy)){
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
results <- list("Hierarchical_levels" = 2,"Order_diversity" = qvalue,"Sample_size" = N, "L1_diversity" = L1_div, "L2_diversity" = L2_div, "Beta_diversity" = beta)
return(results)
}

########### NEED TO BE CORRECTED - 2018/12/16
#Function for 3-level hierarchy
#Samples are weighed to give the same weight to the L2 groups.
###########

#Pre-processing
if(ncol(hierarchy) == 2){
hierarchy[,1] <- as.character(hierarchy[,1])
hierarchy[,2] <- as.character(hierarchy[,2])
if(identical(sort(colnames(otutable)),sort(hierarchy[,1])) == FALSE) stop("OTU names in the OTU table and the hierarchy table do not match")
colnames(hierarchy) <- c("L1","L2")

#L2 data processing
otutable.L2 <- merge(t(otutable),hierarchy, by.x="row.names",by.y="L1")
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- aggregate(otutable.L2[,-ncol(otutable.L2)], by=list(otutable.L2[,ncol(otutable.L2)]), FUN=sum)
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- t(sweep(otutable.L2, 1, rowSums(otutable.L2), FUN="/"))

#Compute L1 weights
weights.L1 <- rep(1,ncol(otutable))/rep(ncol(otutable.L2),ncol(otutable))/rep(table(hierarchy[,2]), table(hierarchy[,2]))

#L1 and L3 diversities
if(missing(tree)){  
  L1_div <- hilldiv::alpha.div(otutable,qvalue,weight=weights.L1)
  L3_div <- hilldiv::gamma.div(otutable,qvalue)
  }else{
  L1_div <- hilldiv::alpha.div(otutable,qvalue,tree,weight=weights.L1)
  L3_div <- hilldiv::gamma.div(otutable,qvalue,tree)
}

#L2 diversity
if(missing(tree)){ 
  L2_div <- hilldiv::alpha.div(otutable.L2,qvalue)
  }else{
  L2_div <- hilldiv::alpha.div(otutable.L2,qvalue,tree)
}

#L1-L2, beta and sample size
beta_1 <- L2_div/L1_div
N1 <- ncol(otutable)

#L2-L3, beta and sample size
beta_2 <- L3_div/L2_div
N2 <- ncol(otutable.L2)

#Return values
results <- list("Hierarchical_levels" = 2,"Order_diversity" = qvalue, "L1_Sample_size" = N1, "L2_Sample_size" = N2, "L1_diversity" = L1_div, "L2_diversity" = L2_div, "L3_diversity" = L3_div, "Beta_diversity_L1_2" = beta_1, "Beta_diversity_L2_3" = beta_2)
return(results)
}

###########
#Function for 4-level hierarchy
###########
  
if(ncol(hierarchy) == 3){
stop("The maximum number of hierarchical levels allowed is 3")
} 
  
#Error if more hierarchical levels
if(ncol(hierarchy) > 3){
stop("The maximum number of hierarchical levels allowed is 3")
} 

}
