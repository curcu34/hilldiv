pair.dis <- function(otutable,qvalue,measure,hierarchy){

#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}
if(missing(measure)) { measure= c("homogeneity","overlap","turnover")}
if(missing(hierarchy)) warning("Assuming no hierarchy")

#Create matrices
L1 <- sort(colnames(otutable))
N <- length(L1)

L1_beta <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_beta) <- L1
rownames(L1_beta) <- L1

if('homogeneity' %in% measure){
L1_homogeneity <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_homogeneity) <- L1
rownames(L1_homogeneity) <- L1
}

if('overlap' %in% measure){
L1_overlap <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_overlap) <- L1
rownames(L1_overlap) <- L1
}

if('turnover' %in% measure){
L1_turnover <- matrix(rep(NA,length(L1)^2), nrow = length(L1), ncol = length(L1))
colnames(L1_turnover) <- L1
rownames(L1_turnover) <- L1
}

for (x in L1){
for (y in L1){
if(is.na(L1_beta[y,x])){ #to avoid repeating mirror operations
combination <- otutable[,c(x,y)]
alpha <- alpha.div(combination,qvalue)
gamma <- gamma.div(combination,qvalue)
beta <- gamma/alpha
L1_beta[x,y] <- beta
results <- list("Beta" = L1_beta)

if('homogeneity' %in% measure){
homogeneity <- ((1/beta) - 1/N)/(1-1/N)
L1_homogeneity[x,y] <- homogeneity
results[["Homogeneity"]] <- L1_homogeneity}

if('overlap' %in% measure){
overlap <-((1/beta)^(1-qvalue) - (1/N)^(1-qvalue)) / (1 - (1/N)^(1-qvalue))
L1_overlap[x,y] <- overlap
results[["Overlap"]] <- L1_overlap}

if('turnover' %in% measure){
turnover <- (beta - 1)/(N-1)
L1_turnover[x,y] <- turnover
results[["Turnover"]] <- L1_turnover}
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
N <- length(L2)

L2_beta <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_beta) <- L2
rownames(L2_beta) <- L2

if('homogeneity' %in% measure){
L2_homogeneity <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_homogeneity) <- L2
rownames(L2_homogeneity) <- L2
}

if('overlap' %in% measure){
L2_overlap <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_overlap) <- L2
rownames(L2_overlap) <- L2
}

if('turnover' %in% measure){
L2_turnover <- matrix(rep(NA,length(L2)^2), nrow = length(L2), ncol = length(L2))
colnames(L2_turnover) <- L2
rownames(L2_turnover) <- L2
}

for (x in L2){
for (y in L2){
if(is.na(L2_beta[y,x])){ #to avoid repeating mirror operations
combination <- otutable.L2[,c(x,y)]
alpha <- alpha.div(combination,qvalue)
gamma <- gamma.div(combination,qvalue)
beta <- gamma/alpha
L2_beta[x,y] <- beta
results[["Beta_L2"]] <- L2_beta

if('homogeneity' %in% measure){
homogeneity <- ((1/beta) - 1/N)/(1-1/N)
L2_homogeneity[x,y] <- homogeneity
results[["Homogeneity_L2"]] <- L2_homogeneity}

if('overlap' %in% measure){
overlap <-((1/beta)^(1-qvalue) - (1/N)^(1-qvalue)) / (1 - (1/N)^(1-qvalue))
L2_overlap[x,y] <- overlap
results[["Overlap_L2"]] <- L2_overlap}

if('turnover' %in% measure){
turnover <- (beta - 1)/(N-1)
L2_turnover[x,y] <- turnover
results[["Turnover_L2"]] <- L2_turnover}
}
}
}
  
}
  
return(results)

}



