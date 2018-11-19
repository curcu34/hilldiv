hill.dis <- function(beta,qvalue,N,metric,type){ 
  
#Quality-check and warnings
if(missing(beta)) stop("Beta diversity value is missing")
if(missing(qvalue)) stop("The order of diversity (q) is missing")
if (qvalue==1) {qvalue=0.99999}
if(missing(N)) stop("The number of samples or groups (N) is missing")
if (missing(metric)) {metric="C"}
if (missing(type)) {type="dissimilarity"}

#Sørensen-type species-overlap (CqN, 1-CqN)
if (metric == "C"){
CqN <- ((1/beta)^(1-qvalue) - (1/N)^(1-qvalue)) / (1 - (1/N)^(1-qvalue))
return(CqN)
if (type == "dissimilarity"){
1_CqN <- 1 - CqN
return(1_CqN)
}
}

#Jaccard-type species-overlap (UqN, 1-UqN)
if (metric == "U"){
if (type == "dissimilarity"){ }
}

#Sørensen-type turnover (VqN, 1-VqN)
if (metric == "V"){
VqN <-  (beta - 1)/(N-1)
return(VqN)
if (type == "similarity"){ 
1_VqN <- 1 - VqN
return(1_VqN)
}
} 
  
#Jaccard-type turnover (SqN, 1-SqN)
if (metric == "S"){
if (type == "dissimilarity"){ }
} 
  
}
