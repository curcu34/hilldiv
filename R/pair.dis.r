pair.dis <- function(otutable,qvalue,tree,weight,hierarchy,level,metric){

#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}\
if(missing(weight)) { weight= rep(1/ncol(otutable),ncol(otutable))}
if(length(weight) != ncol(otutable)) stop("The length of the weight vector does not match the number of sequences")
names(weight) = colnames(otutable)
if(!missing(hierarchy)) {
if(missing(level)) {level= c(1:ncol(hierarchy))}
}
if(missing(hierarchy)) {level=1}
if(missing(metric)) { metric= c("C","U","V","S")}

#Declare fast alpha and gamma phylodiversities (without ltips, as it is the same for all combinations)

alpha.div.fast <- function(otutable,qvalue,tree,weight){
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
if(ape::is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")
if (qvalue==1) {qvalue=0.99999}
otutable <- as.data.frame(otutable)
wj <- weight
N <- ncol(otutable)
Li <- tree$edge.length
aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(otutable[TipVector,]))), ncol = N, byrow = TRUE)
aij.wj <- sweep(aij, 2, wj, "*")
T <- sum(sweep(aij.wj, 1, Li, "*"))
L <- matrix(rep(Li, N), ncol = N)
wm <-  matrix(rep(wj, length(Li)), ncol = N,byrow=TRUE)
i <-  which(aij > 0)
phylodiv <- sum(L[i] * (aij[i]*wm[i]/T)^qvalue)^(1/(1 - qvalue))/(N*T)
return(phylodiv)
}

gamma.div.fast <- function(otutable,qvalue,tree,weight){
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999} # change q to the limit of the unity (0.99999) if q=1
if(ape::is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")
if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
otutable <- as.data.frame(otutable)
wj <- weight
N <- ncol(otutable)
Li <- tree$edge.length
aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(otutable[TipVector,]))), ncol = N, byrow = TRUE)
aij.wj <- sweep(aij, 2, wj, "*")
ai <- rowSums(aij.wj)
T <- sum(sweep(aij.wj, 1, Li, "*"))
L <- matrix(rep(Li, N), ncol = N)
Li <- Li[ai != 0] #Remove zeros
ai <- ai[ai != 0] #Remove zeros
wm <-  matrix(rep(wj, length(Li)), ncol = N, byrow=TRUE)
phylodiv <- (sum(Li * (ai/T)^qvalue)^(1/(1 - qvalue)))/T
return(phylodiv)
}

#Generate ltips
if(!missing(tree)){
ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
}

#####
#  First herarchical level
#####

if("1" %in% level){

#Create matrices
L1 <- sort(colnames(otutable))

L1_beta <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_beta) <- L1
rownames(L1_beta) <- L1

if('C' %in% metric){
L1_CqN <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_CqN) <- L1
rownames(L1_CqN) <- L1
}

if('U' %in% metric){
L1_UqN <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_UqN) <- L1
rownames(L1_UqN) <- L1
}

if('V' %in% metric){
L1_VqN <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_VqN) <- L1
rownames(L1_VqN) <- L1
}

if('S' %in% metric){
L1_SqN <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_SqN) <- L1
rownames(L1_SqN) <- L1
}

#Fill matrices

for (x in L1){
for (y in L1){
if(is.na(L1_beta[x,y])){ #to avoid repeating mirror operations
combination <- otutable[,c(y,x)]
combination.weight <- weight[c(y,x)]
combination.weight <- combination.weight/(sum(combination.weight))

if(identical(combination[,1],combination[,2]) == TRUE){
    beta <- NA
}else{
    if(missing(tree)){
    alpha <- hilldiv::alpha.div(combination,qvalue,weight=combination.weight)
    gamma <- hilldiv::gamma.div(combination,qvalue,weight=combination.weight)
    }else{
    alpha <- alpha.div.fast(combination,qvalue,tree,weight=combination.weight)
    gamma <- gamma.div.fast(combination,qvalue,tree,weight=combination.weight)
    }
    beta <- gamma/alpha
}

L1_beta[y,x] <- beta
results <- list("L1_Beta" = L1_beta)

if('C' %in% metric){
disC <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=2,metric="C",type="dissimilarity")
L1_CqN[y,x] <- disC
results[["L1_CqN"]] <- L1_CqN}

if('U' %in% metric){
disU <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=2,metric="U",type="dissimilarity")
L1_UqN[y,x] <- disU
results[["L1_UqN"]] <- L1_UqN}

if('V' %in% metric){
disV <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=2,metric="V",type="dissimilarity")
L1_VqN[y,x] <- disV
results[["L1_VqN"]] <- L1_VqN}

if('S' %in% metric){
disS <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=2,metric="S",type="dissimilarity")
L1_SqN[y,x] <- disS
results[["L1_SqN"]] <- L1_SqN}

}
}
}
}

#If hierarchy specified, also do:
if("2" %in% level){

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

#Prepare L2 weights
weight.L2pre <- merge(t(t(weight)),hierarchy, by.x="row.names",by.y="L1")
weight.L2pre <- aggregate(weight.L2[,2], by=list(weight.L2[,3]), FUN=sum)
weight.L2 <- weight.L2pre[,2]
names(weight.L2) <- weight.L2pre[,1]

#Create L2 matrices
L2 <- sort(colnames(otutable.L2))

L2_beta <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_beta) <- L2
rownames(L2_beta) <- L2

if('C' %in% metric){
L2_CqN <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_CqN) <- L2
rownames(L2_CqN) <- L2
}

if('U' %in% metric){
L2_UqN <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_UqN) <- L2
rownames(L2_UqN) <- L2
}

if('V' %in% metric){
L2_VqN <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_VqN) <- L2
rownames(L2_VqN) <- L2
}

if('S' %in% metric){
L2_SqN <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_SqN) <- L2
rownames(L2_SqN) <- L2
}

#Fill L2 matrices

for (x in L2){
for (y in L2){
if(is.na(L2_beta[x,y])){ #to avoid repeating mirror operations
combination <- otutable.L2[,c(y,x)]
combination.weight <- weight.L2[c(y,x)]
combination.weight <- combination.weight/(sum(combination.weight))

if(identical(combination[,1],combination[,2]) == TRUE){
    beta <- NA
}else{
    if(missing(tree)){
    alpha <- hilldiv::alpha.div(combination,qvalue,weight=combination.weight)
    gamma <- hilldiv::gamma.div(combination,qvalue,weight=combination.weight)
    }else{
    alpha <- alpha.div.fast(combination,qvalue,tree,weight=combination.weight)
    gamma <- gamma.div.fast(combination,qvalue,tree,weight=combination.weight)
    }
    beta <- gamma/alpha
}

L2_beta[y,x] <- beta
if(exists("results")){
results[["L2_Beta"]] <- L2_beta
}else{
results <- list("L2_Beta" = L2_beta)
}

if('C' %in% metric){
disC <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=2,metric="C",type="dissimilarity")
L2_CqN[y,x] <- disC
results[["L2_CqN"]] <- L2_CqN}

if('U' %in% metric){
disU <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=2,metric="U",type="dissimilarity")
L2_UqN[y,x] <- disU
results[["L2_UqN"]] <- L2_UqN}

if('V' %in% metric){
disV <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=2,metric="V",type="dissimilarity")
L2_VqN[y,x] <- disV
results[["L2_VqN"]] <- L2_VqN}

if('S' %in% metric){
disS <- hilldiv::beta.dis(beta=beta,qvalue=qvalue,N=2,metric="S",type="dissimilarity")
L2_SqN[y,x] <- disS
results[["L2_SqN"]] <- L2_SqN}

}
}
}

}

return(results)

}
