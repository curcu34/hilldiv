beta.dis <- function(beta,qvalue,N,metric,type){ 
  
#Quality-check and warnings
if(missing(beta)) stop("Beta diversity value is missing")
if(missing(qvalue)) stop("The order of diversity (q) is missing")
if (qvalue==1) {qvalue=0.99999}
if(missing(N)) stop("The number of samples or groups (N) is missing")
if (missing(metric)) {metric="C"}
if (missing(type)) {type="dissimilarity"}

#Sørensen-type overlap (CqN, 1-CqN)
if (metric == "C"){
CqN <- ((1/beta)^(qvalue-1) - (1/N)^(qvalue-1)) / (1 - (1/N)^(qvalue-1))
  if (type == "dissimilarity"){
  rCqN <- 1 - CqN
  return(rCqN)
  }
return(CqN)
}

#Jaccard-type overlap (UqN, 1-UqN)
if (metric == "U"){
UqN <- ((1/beta)^(1-qvalue) - (1/N)^(1-qvalue)) / (1 - (1/N)^(1-qvalue))
  if (type == "dissimilarity"){
  rUqN <- 1 - UqN
  return(rUqN)
  }
return(UqN)
}

#Sørensen-type turnover (VqN, 1-VqN)
if (metric == "V"){
VqN <-  (beta - 1)/(N-1)
  if (type == "similarity"){ 
  rVqN <- 1 - VqN
  return(rVqN)
  }
return(VqN)
} 
  
#Jaccard-type turnover (SqN, 1-SqN)
if (metric == "S"){
SqN <- ((1/beta) - 1/N)/(1-1/N)
  if (type == "dissimilarity"){
  rSqN <- 1 - SqN
  return(rSqN)
  }
return(SqN)
} 
  
}
