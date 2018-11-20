pair.dis <- function(otutable,qvalue,tree,hierarchy,measure){

#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}
if(missing(measure)) { measure= c("C","U","V","S")}
if(missing(hierarchy)) warning("Assuming no hierarchy")

#Create matrices
L1 <- sort(colnames(otutable))
N <- 2
  
L1_beta <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_beta) <- L1
rownames(L1_beta) <- L1

if('C' %in% measure){
L1_CqN <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_CqN) <- L1
rownames(L1_CqN) <- L1
}

if('U' %in% measure){
L1_UqN <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_UqN) <- L1
rownames(L1_UqN) <- L1
}

if('V' %in% measure){
L1_VqN <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_VqN) <- L1
rownames(L1_VqN) <- L1
}

if('S' %in% measure){
L1_SqN <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_SqN) <- L1
rownames(L1_SqN) <- L1
}

#Fill matrices  
  
for (x in L1){
for (y in L1){
if(is.na(L1_beta[x,y])){ #to avoid repeating mirror operations
combination <- otutable[,c(y,x)]

if(identical(combination[,1],combination[,2]) == TRUE){
    beta <- NA
}else{
    if(missing(tree)){
    alpha <- hilldiv::alpha.div(combination,qvalue)
    gamma <- hilldiv::gamma.div(combination,qvalue)
    }else{
    alpha <- hilldiv::alpha.div(combination,qvalue,tree)
    gamma <- hilldiv::gamma.div(combination,qvalue,tree)
    }
    beta <- gamma/alpha
}
  
L1_beta[y,x] <- beta
results <- list("Beta" = L1_beta)
  
if('C' %in% measure){
disC <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="C",type="dissimilarity")
L1_CqN[y,x] <- disC
results[["1_CqN"]] <- L1_CqN}
  
if('U' %in% measure){
disU <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="U",type="dissimilarity")
L1_UqN[y,x] <- disU
results[["1_UqN"]] <- L1_UqN}

if('V' %in% measure){
disV <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="V",type="dissimilarity")
L1_VqN[y,x] <- disV
results[["1_VqN"]] <- L1_VqN}

if('S' %in% measure){
disS <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="S",type="dissimilarity")
L1_SqN[y,x] <- disS
results[["1_SqN"]] <- L1_SqN}
  
}
}
}

#If hierarchy specified, also do:
if(ncol(hierarchy) == 2){

#Check hierarchy file
hierarchy[,1] <- as.character(hierarchy[,1])
hierarchy[,2] <- as.character(hierarchy[,2])
if(identical(sort(colnames(otutable)),sort(hierarchy[,1])) == FALSE) stop("OTU names in the OTU table and the hierarchy table do not match")
colnames(hierarchy) <- c("L1","L2")

#Prepare L2 OTU table
otutable.L2 <- merge(t(otutable),hierarchy, by.x="row.names",by.y="L1")
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- aggregate(otutable.L2[,-ncol(otutable.L2)], by=list(otutable.L2[,ncol(otutable.L2)]), FUN=sum)
rownames(otutable.L2) <- otutable.L2[,1]
otutable.L2 <- otutable.L2[,-1]
otutable.L2 <- t(sweep(otutable.L2, 1, rowSums(otutable.L2), FUN="/"))

#Create L2 matrices
L2 <- sort(colnames(otutable.L2))
N <- 2

L2_beta <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_beta) <- L2
rownames(L2_beta) <- L2

if('C' %in% measure){
L2_CqN <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_CqN) <- L2
rownames(L2_CqN) <- L2
}

if('U' %in% measure){
L2_UqN <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_UqN) <- L2
rownames(L2_UqN) <- L2
}

if('V' %in% measure){
L2_VqN <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_VqN) <- L2
rownames(L2_VqN) <- L2
}

if('S' %in% measure){
L2_SqN <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_SqN) <- L2
rownames(L2_SqN) <- L2
}
  
#Fill L2 matrices

for (x in L2){
for (y in L2){
if(is.na(L2_beta[x,y])){ #to avoid repeating mirror operations
combination <- otutable.L2[,c(y,x)]
  
if(identical(combination[,1],combination[,2]) == TRUE){
    beta <- NA
}else{
    if(missing(tree)){
    alpha <- hilldiv::alpha.div(combination,qvalue)
    gamma <- hilldiv::gamma.div(combination,qvalue)
    }else{
    alpha <- hilldiv::alpha.div(combination,qvalue,tree)
    gamma <- hilldiv::gamma.div(combination,qvalue,tree)
    }
    beta <- gamma/alpha
}
  
L2_beta[y,x] <- beta
results[["Beta_L2"]] <- L2_beta

if('C' %in% measure){
disC <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="C",type="dissimilarity")
L2_CqN[y,x] <- disC
results[["1_CqN"]] <- L2_CqN}
  
if('U' %in% measure){
disU <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="U",type="dissimilarity")
L2_UqN[y,x] <- disU
results[["1_UqN"]] <- L2_UqN}

if('V' %in% measure){
disV <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="V",type="dissimilarity")
L2_VqN[y,x] <- disV
results[["1_VqN"]] <- L2_VqN}

if('S' %in% measure){
disS <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=N,metric="S",type="dissimilarity")
L2_SqN[y,x] <- disS
results[["1_SqN"]] <- L2_SqN}
 
}
}
}
  
return(results)

}
