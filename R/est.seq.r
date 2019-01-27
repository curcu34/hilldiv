est.seq <- function(otutable,qvalue,tree,colour,size){
  
#WARNINGS & DEFAULTS
if((qvalue != 0) & (qvalue != 1) & (qvalue != 2))  stop("The order of diversity (q) must to be 0, 1 or 2.")
maxsize=max(colSums(otutable))  
if(missing(size)){size=seq(1,maxsize*2,round(maxsize*2/20))}
  
#SET ANALYSIS TYPE
if(missing(tree)){

est.inext <- iNEXT(otutable[,1], q=qvalue, datatype="abundance")

}else{

}
  
}
