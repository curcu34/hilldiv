beta.dis <- function(beta,qvalue,N,metric,type){ 
  
#Identify input type  
if(missing(beta)) stop("Beta diversity value or div.part output object is missing")
if(class(beta) == "numeric"){input="beta"}
if(class(beta) == "list"){
input="list"
if(names(beta[1]) != "Order_diversity") stop("The input object is not valid")
qvalue <- beta[[2]]
N <- beta[[3]]
beta <- beta[[6]]
}
  
#Quality-check and warnings
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

#Sørensen-type turnover-complement (VqN, 1-VqN)
if (metric == "V"){
VqN <-  (N - beta)/(N-1)
  if (type == "dissimilarity"){ 
  rVqN <- 1 - VqN
  return(rVqN)
  }
return(VqN)
} 
  
#Jaccard-type turnover-complement (SqN, 1-SqN)
if (metric == "S"){
SqN <- ((1/beta) - 1/N)/(1-1/N)
  if (type == "dissimilarity"){
  rSqN <- 1 - SqN
  return(rSqN)
  }
return(SqN)
} 
  
}
