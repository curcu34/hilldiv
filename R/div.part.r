div.part <- function(otutable,qvalue,hierarchy,tree) {
  
#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}
if(missing(hierarchy)) warning("Assuming a two-level hierarchy: 1) sample, 2) total dataset")

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

#Dissimilarity measures
disC <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="C",type="dissimilarity")
disU <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="U",type="dissimilarity")
disV <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="V",type="dissimilarity")
disS <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="S",type="dissimilarity")

#Return values
results <- list("Sample_size" = N, "L1_diversity" = L1_div, "L2_diversity" = L2_div, "Beta_diversity" = beta, "1_CqN" = disC, "1_UqN" = disU, "1_VqN" = disV, "1_SqN" = disS)
return(results)
}

###########
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
  
#L1-L2, Dissimilarities
disC_1 <- hilldiv::beta.dis(beta=beta_1,qvalue=qvalue,N=N1,metric="C",type="dissimilarity")
disU_1 <- hilldiv::beta.dis(beta=beta_1,qvalue=qvalue,N=N1,metric="U",type="dissimilarity")
disV_1 <- hilldiv::beta.dis(beta=beta_1,qvalue=qvalue,N=N1,metric="V",type="dissimilarity")
disS_1 <- hilldiv::beta.dis(beta=beta_1,qvalue=qvalue,N=N1,metric="S",type="dissimilarity")

#L2-L3, beta and sample size
beta_2 <- L3_div/L2_div
N2 <- ncol(otutable.L2)

#L2-L3, Dissimilarities
disC_2 <- hilldiv::beta.dis(beta=beta_2,qvalue=qvalue,N=N2,metric="C",type="dissimilarity")
disU_2 <- hilldiv::beta.dis(beta=beta_2,qvalue=qvalue,N=N2,metric="U",type="dissimilarity")
disV_2 <- hilldiv::beta.dis(beta=beta_2,qvalue=qvalue,N=N2,metric="V",type="dissimilarity")
disS_2 <- hilldiv::beta.dis(beta=beta_2,qvalue=qvalue,N=N2,metric="S",type="dissimilarity")

#Return values
results <- list("L1_size" = N1, "L2_size" = N2, "L1_diversity" = L1_div, "L2_diversity" = L2_div, "L3_diversity" = L3_div, "Beta_diversity_L1_2" = beta_1, "Beta_diversity_L2_3" = beta_2,  "1_CqN_L1_L2" = disC_1, "1_UqN_L1_L2" = disU_1,"VqN_L1_L2" = disV_1, "1_SqN_L1_L2" = disS_1, "1_CqN_L2_L3" = disC_2, "1_UqN_L2_L3" = disU_2,"VqN_L2_L3" = disV_2, "1_SqN_L2_L3" = disS_2)
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
