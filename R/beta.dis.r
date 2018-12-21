beta.dis <- function(beta,qvalue,N,metric,type){ 
  
#Identify input type  
if(missing(beta)) stop("Beta diversity value or div.part output object is missing")
if(class(beta) == "numeric"){input="beta"}
if(class(beta) == "list"){
input="list"
if(names(beta[3]) != "Order_diversity") stop("The input object is not valid")
  if(beta[[1]] == 2){
    qvalue <- beta[[3]]
    N <- beta[[4]]
    beta <- beta[[7]]
  }
  if(beta[[1]] == 3){
    qvalue <- beta[[3]]
    N1 <- beta[[4]]
    N2 <- beta[[5]]
    beta1 <- beta[[9]]
    beta2 <- beta[[10]]
  }
}
  
#Quality-check and warnings
if(missing(qvalue)) stop("The order of diversity (q) is missing")
if (qvalue==1) {qvalue=0.99999}
if(missing(N) && missing(N1)) stop("The number of samples or groups (N) is missing")
if (missing(metric)) {metric="U"}
if (missing(type)) {type="dissimilarity"}

###### MULTIPLE HIERARCHIES NEED TO BE ADDED, and similarity functions updated - 2018/12/21  
  
#Sørensen-type overlap (CqN, 1-CqN)
if (metric == "C"){
  if(exists("N1")){
  CqN1 <- CqN(beta1,qvalue,N1)
  CqN2 <- CqN(beta2,qvalue,N2)
  mCqN <- list(CqN1,CqN2)
  if (type == "dissimilarity"){
    rCqN1 <- 1 - CqN1
    rCqN2 <- 1 - CqN2
    mrCqN <- list(rCqN1,rCqN2)
    return(mrCqN)
  }
  return(mCqN)  
  }else{
  CqN <- CqN(beta,qvalue,N)
  if (type == "dissimilarity"){
  rCqN <- 1 - CqN
  return(rCqN)
  }
  return(CqN)
}
  
}

#Jaccard-type overlap (UqN, 1-UqN)
if (metric == "U"){
  if(exists("N1")){
  UqN1 <- UqN(beta1,qvalue,N1)
  UqN2 <- UqN(beta2,qvalue,N2)
  mUqN <- list(UqN1,UqN2)
  if (type == "dissimilarity"){
    rUqN1 <- 1 - UqN1
    rUqN2 <- 1 - UqN2
    mrUqN <- list(rUqN1,rUqN2)
    return(mrUqN)
  }
  return(mUqN)  
  }else{
  UqN <- UqN(beta,qvalue,N)
  if (type == "dissimilarity"){
  rUqN <- 1 - UqN
  return(rUqN)
  }
  return(UqN)
}
  
}

#Sørensen-type turnover-complement (VqN, 1-VqN)
if (metric == "V"){
  if(exists("N1")){
  VqN1 <- VqN(beta1,N1)
  VqN2 <- VqN(beta2,N2)
  mVqN <- list(VqN1,VqN2)
  if (type == "dissimilarity"){
    rVqN1 <- 1 - VqN1
    rVqN2 <- 1 - VqN2
    mrVqN <- list(rVqN1,rVqN2)
    return(mrVqN)
  }
  return(mVqN)  
  }else{
  VqN <- VqN(beta,N)
  if (type == "dissimilarity"){
  rVqN <- 1 - VqN
  return(rVqN)
  }
  return(VqN)
}

} 
  
#Jaccard-type turnover-complement (SqN, 1-SqN)
if (metric == "S"){
  if(exists("N1")){
  SqN1 <- SqN(beta1,N1)
  SqN2 <- SqN(beta2,N2)
  mSqN <- list(SqN1,SqN2)
  if (type == "dissimilarity"){
    rSqN1 <- 1 - SqN1
    rSqN2 <- 1 - SqN2
    mrSqN <- list(rSqN1,rSqN2)
    return(mrSqN)
  }
  return(mSqN)  
  }else{
  SqN <- SqN(beta,N)
  if (type == "dissimilarity"){
  rSqN <- 1 - SqN
  return(rSqN)
  }
  return(SqN)
}
  
}
