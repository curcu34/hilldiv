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
if (type == "dissimilarity"){ }
}

#Jaccard-type species-overlap (UqN, 1-UqN)
if (metric == "U"){
if (type == "dissimilarity"){ }
}

#Sørensen-type turnover (VqN, 1-VqN)
if (metric == "V"){
if (type == "dissimilarity"){ }
} 
  
#Jaccard-type turnover (SqN, 1-SqN)
if (metric == "S"){
if (type == "dissimilarity"){ }
} 
  
}
