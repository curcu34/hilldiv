pair.dis <- function(otutable,qvalue,hierarchy,measure){

#Quality-check and warnings
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}
if(missing(hierarchy)) warning("Assuming no hierarchy")
if(missing(measure)) { measure= c("homogeneity","overlap","turnover")}

##Function for no hierarchy
if(missing(hierarchy)){

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

return(results)
}

#Function for 2-level hierarchy

  
}



