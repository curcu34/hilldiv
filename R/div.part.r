div.part <- function(otutable,qvalue,hierarchy,tree,type) {
  
#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(qvalue==1) {qvalue=0.99999}
if(missing(type)){type == "abundance"}
if(type == "incidence"){
  if(missing(hierarchy)) stop("Diversity partitioning based on incidence data requires a hierarchy table")
}  
###########
#Function for 2-level hierarchy
########### 

#L1 and L2 diversities
if(missing(hierarchy)){
if(type == "incidence") stop("Diversity partitioning based on incidence data requires a hierarchy table")
if(type == "estimate") stop("Estimated diversity partitioning based on incidence data requires a hierarchy table")
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
                           
}else if(type == "estimate"){
#####  
# Incidence-estimation
#####
  
if((qvalue != 0) & (qvalue != 1) & (qvalue != 2))  stop("For estimated diversity partitioning the order of diversity (q) must to be 0, 1 or 2.")

#Transform to iNEXT/iNextPD format  
otutable.inext.L2 <- to.inext(otutable,hierarchy=sample.species,type="incidence_raw")
otutable.inext.L3 <- to.inext(otutable,type="incidence_raw")

if(missing(tree)){  
  #Alpha diversity
  if(qvalue == 0){L2_div <- mean(ChaoSpecies(otutable.inext.L2,datatype="incidence_raw", conf=0.95)[,2])}
  if(qvalue == 1){L2_div <- exp(mean(ChaoEntropy(otutable.inext.L2,datatype="incidence_raw", transform=FALSE, conf=0.95)[,2]))}
  if(qvalue == 2){L2_div <- 1/(1-mean(EstSimpson(otutable.inext.L2,datatype="incidence_raw", transform=FALSE, conf=0.95)[,2]))}
  #Gamma diversity  
  if(qvalue == 0){L3_div <- ChaoSpecies(otutable.inext.L3,datatype="incidence_raw", conf=0.95)[,2]}
  if(qvalue == 1){L3_div <- ChaoEntropy(otutable.inext.L3,datatype="incidence_raw", transform=TRUE, conf=0.95)[,2]}
  if(qvalue == 2){L3_div <- EstSimpson(otutable.inext.L3,datatype="incidence_raw", transform=TRUE, conf=0.95)[,2]}

}else{
  tree.phylog <- phylo.to.phylog(tree)
  #Alpha diversity
  if(qvalue == 0){L2_div <- mean(estPD(otutable.inext.L2,labels=rownames(otutable),phy=tree.phylog,q=0,datatype="incidence_raw", se=FALSE, conf=0.95)[,2])}
  if(qvalue == 1){L2_div <- exp(mean(log(estPD(otutable.inext.L2,labels=rownames(otutable),phy=tree.phylog,q=1,datatype="incidence_raw", se=FALSE, conf=0.95)[,2])))}
  if(qvalue == 2){
    L2_div.raw <- estPD(otutable.inext.L2,labels=rownames(otutable),phy=tree.phylog,q=2,datatype="incidence_raw", se=FALSE, conf=0.95)[,2]}
    L2_div.raw2 <- mean(-(1-L2_div.raw)/L2_div.raw)
    L2_div <- 1/(1-L2_div.raw2)}
  #Gamma diversity  
    L3_div <- estPD(otutable.inext.L3,labels=rownames(otutable),phy=tree.phylog,q=qvalue,datatype="incidence_raw", se=FALSE, conf=0.95)
  }
  #L2-L3, beta and sample size
  beta <- L3_div/L2_div
  N <- length(otutable.inext.L2)
  
  #Return values
  results <- list("Hierarchical_levels" = 2,"Type" = type,"Order_diversity" = qvalue,"Sample_size" = N, "L2_diversity" = L2_div, "L3_diversity" = L3_div, "Beta_diversity" = beta)
  return(results)
                                             
}else{
stop("The type of diversity partition provided is incorrect.")
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
