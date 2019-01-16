gamma.est <- function(otutable,qvalue,hierarchy,tree,steps,iter,pred.size,summary,conf){

if(missing(hierarchy) == TRUE){
##############################
# HIERARCHY IS NOT SPECIFIED #
##############################

cat("Rarefying gamma diversities \n")

max.steps <- ncol(otutable)
selected.steps <- max.steps * steps / 100
interval <- max.steps/selected.steps
step.vector <- seq(from = 1, to = max.steps, by = interval)
if(step.vector[length(step.vector)] != max.steps){ step.vector <- c(step.vector,max.steps)} #Add last number if is not included	

## Rarefy gamma diversities 
matrix <- c()
for (i in 1:iter){#open iteration
  vector <- c()
  for (s in step.vector){#open steps
      subset <- sample(max.steps, s, replace = FALSE)
      otutable.subset <- otutable[,subset]
      if(is.null(dim(otutable.subset)) == TRUE){#If 1 sample, calculate Hill number
	    	names(otutable.subset) <- rownames(otutable)
          value <- hilldiv::hill.div(otutable.subset,qvalue=qvalue)
	  	}else{#If >1 sample, calculate Gamma diversity
          value <- hilldiv::gamma.div(otutable.subset,qvalue=qvalue)
	  	}
    vector <- c(vector,value)
    }#close steps
matrix <- rbind(matrix,vector)
}#close iteration

#Melt matrix (MODIFY!!!!!)
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]
colnames(matrix) <- step.vector
matrix.melted <- melt(matrix)
colnames(matrix.melted) <- c("Step","Value")
matrix.melted[,c(1,2)] <- apply(matrix.melted[,c(1,2)], 2, function(x) as.numeric(as.character(x))) #change to numeric

#Observed
raref.obs <- aggregate(matrix.melted[,2],by=list(matrix.melted[,1]),FUN=mean)
raref.error <- aggregate(matrix.melted[,2],by=list(matrix.melted[,1]),function(x) qnorm(0.975)*sd(x)/sqrt(length(x)))
raref.lwr <- raref.obs[,2] - raref.error[,2]
raref.upr <- raref.obs[,2] + raref.error[,2]
raref.obs <- cbind(raref.obs,raref.lwr,raref.upr)
colnames(raref.obs) <- c("Step","Observed","Obs.lwr","Obs.upr")

#Generate NA-s
pred.steps <- seq((nrow(raref.obs)+1),pred.size)
raref.NA <- cbind(pred.steps,matrix(rep(NA,length(pred.steps)*3), ncol = 3))
colnames(raref.NA) <- c("Step","Observed","Obs.lwr","Obs.upr")
raref.obs.NA <- rbind(raref.obs,raref.NA)

#Model
mod <- nls(Value ~ SSasymp(Step, Asym, resp0, lrc), data = matrix.melted)

#Asymptote
max.obs <- tail(matrix.melted$Value,1)
asymp.est <- mod$m$getPars()[1]
coverage <- max.obs/asymp.est
correlation <- cor(matrix.group$Value,predict(mod))

#Predict
preds = data.frame(Step = sort(rep(c(1:pred.size), times = iter)))
preds$Model = predict(mod, newdata = preds)

#Calculate necessary sample size and completeness
Modelled <- unique(preds$Model)
if(max(Modelled) > (asymp.est*0.99)){
est.sample.size <- min(which(Modelled > (asymp.est*0.99)))
completeness <- nrow(raref.obs)/est.sample.size
if(completeness > 1){completeness = 1}
}else{
est.sample.size <- "The sample size specified for gamma diversity prediction was too low to estimate a sample size. Choose a larger value."
completeness <- "The sample size specified for gamma diversity prediction was too low to estimate sampling completeness. Choose a larger value."
}

#Create coverage and completeness vectors
coverage.vector <- round(Modelled/(asymp.est*0.99),4)
coverage.vector[which(coverage.vector > 1)] <- 1
if(max(Modelled) > (asymp.est*0.99)){
completeness.vector <- round(raref.obs.NA$Step/est.sample.size,4)
completeness.vector[which(completeness.vector > 1)] <- 1
}else{
completeness.vector <- rep(NA,length(coverage.vector))
}

if(conf == TRUE){
#Bootstrap
set.seed(42)
myboot <- boot(matrix.group, function(matrix.group, ind, mod, preds) {
  matrix.group$Value <- fitted(mod) + residuals(mod)[ind]
  tryCatch(predict(nls(Value ~ SSasymp(Step, Asym, resp0, lrc), data = matrix.group), newdata = preds), error = function(e) preds$Step * NA)}, mod = mod, preds = preds, R = 1e4)

#Calculate confidence intervals
CI <- t(sapply(seq_len(nrow(preds)), function(i) boot.ci(myboot, type = "bca", index = i)$bca[4:5]))
colnames(CI) <- c("Mod.lwr","Mod.upr")
result <- cbind(raref.obs.NA,Modelled, unique(CI))
}else{
result <- cbind(raref.obs.NA[,c(1:2)],Modelled,coverage.vector,completeness.vector)
colnames(group.result)[4] <- c("Coverage")
colnames(group.result)[5] <- c("Completeness")

}

results[["Steps"]] <- result
results[["Model.fit"]] <- correlation
results[["Asymptote"]] <- unname(asymp.est)
results[["Coverage"]] <- unname(coverage)
results[["Min.sample"]] <- est.sample.size
results[["Completeness"]] <- completeness
return(results)

}else{
##########################
# HIERARCHY IS SPECIFIED #
##########################

max.steps <- max(table(hierarchy[,2]))
selected.steps <- max.steps * steps / 100
interval <- max.steps/selected.steps
step.vector <- seq(from = 1, to = max.steps, by = interval)
if(step.vector[length(step.vector)] != max.steps){ step.vector <- c(step.vector,max.steps)} #Add last number if is not included	

#Define groups	
groups <- as.character(unique(hierarchy[,2]))

## Rarefy gamma diversities 
if(missing(tree)){ 
cat("Rarefying gamma diversities \n")
}else{
cat("Rarefying gamma phylodiversities\nNOTE: this step might take hours in the case of large phylogenetic trees \n")
}
#Loop across steps, groups and number of iterations	
matrix <- c()
for (g in groups){#open group
cat("    ",g,"\n")

samples <- as.character(hierarchy[which(hierarchy[,2] == g),1])
otutable.group <- otutable[,samples]
    for (i in 1:iter){#open iteration
	  vector <- c()
	  for (s in step.vector){#open steps
		if(s > ncol(otutable.group)){#if sample size is smaller  
		value <- NA}else{
	      subset <- sample(ncol(otutable.group), s, replace = FALSE)
	      otutable.subset <- otutable.group[,subset]
	      if(is.null(dim(otutable.subset)) == TRUE){#If 1 sample, calculate Hill number
		  names(otutable.subset) <- rownames(otutable)
			  if(missing(tree)){ 
			  value <- hilldiv::hill.div(otutable.subset,qvalue=qvalue)
			  }else{
			  value <- hilldiv::hill.div(otutable.subset,qvalue=qvalue,tree=tree)
			  }
		  }else{#If >1 sample, calculate Gamma diversity
			if(missing(tree)){ 
			  value <- hilldiv::gamma.div(otutable.subset,qvalue=qvalue)
			  }else{
			  value <- hilldiv::gamma.div(otutable.subset,qvalue=qvalue,tree=tree)
			  }
		      }
		}
	    vector <- c(vector,value)
	    }#close steps
	 vector <- c(g,vector)
	 matrix <- rbind(matrix,vector)    
    }#close iteration
}#close group  

#Melt matrix
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]
colnames(matrix) <- step.vector
matrix.melted <- melt(matrix)
colnames(matrix.melted) <- c("Group","Step","Value")

## Obtain asymptote estimations and functions

cat("Performing diversity estimations \n")
all.results <- c()

for (g in groups){#open group
if(conf == TRUE){
cat("    ",g,"\n")
}

matrix.group <- matrix.melted[which(matrix.melted$Group == g),c(2,3)] #subset group
matrix.group <- matrix.group[rowSums(is.na(matrix.group)) <= 0,] # remove NAs
matrix.group[,c(1,2)] <- apply(matrix.group[,c(1,2)], 2, function(x) as.numeric(as.character(x))) #change to numeric

#Observed
raref.obs <- aggregate(matrix.group[,2],by=list(matrix.group[,1]),FUN=mean)
raref.error <- aggregate(matrix.group[,2],by=list(matrix.group[,1]),function(x) qnorm(0.975)*sd(x)/sqrt(length(x)))
raref.lwr <- raref.obs[,2] - raref.error[,2]
raref.upr <- raref.obs[,2] + raref.error[,2]
raref.obs <- cbind(raref.obs,raref.lwr,raref.upr)
colnames(raref.obs) <- c("Step","Observed","Obs.lwr","Obs.upr")

#Generate NA-s
pred.steps <- seq((nrow(raref.obs)+1),pred.size)
raref.NA <- cbind(pred.steps,matrix(rep(NA,length(pred.steps)*3), ncol = 3))
colnames(raref.NA) <- c("Step","Observed","Obs.lwr","Obs.upr")
raref.obs.NA <- rbind(raref.obs,raref.NA)

#Model
mod <- nls(Value ~ SSasymp(Step, Asym, resp0, lrc), data = matrix.group)

#Asymptote
max.obs <- tail(matrix.group$Value,1)
asymp.est <- mod$m$getPars()[1]
coverage <- max.obs/asymp.est
correlation <- cor(matrix.group$Value,predict(mod))

#Predict
preds = data.frame(Step = sort(rep(c(1:pred.size), times = iter)))
preds$Model = predict(mod, newdata = preds)

#Calculate necessary sample size and completeness
Modelled <- unique(preds$Model)
if(max(Modelled) > (asymp.est*0.99)){
est.sample.size <- min(which(Modelled > (asymp.est*0.99)))
completeness <- nrow(raref.obs)/est.sample.size
if(completeness > 1){completeness = 1}
}else{
est.sample.size <- "The sample size specified for gamma diversity prediction was too low to estimate a sample size. Choose a larger value."
completeness <- "The sample size specified for gamma diversity prediction was too low to estimate sampling completeness. Choose a larger value."
}

#Create coverage and completeness vectors
coverage.vector <- round(Modelled/(asymp.est*0.99),4)
coverage.vector[which(coverage.vector > 1)] <- 1
if(max(Modelled) > (asymp.est*0.99)){
completeness.vector <- round(raref.obs.NA$Step/est.sample.size,4)
completeness.vector[which(completeness.vector > 1)] <- 1
}else{
completeness.vector <- rep(NA,length(coverage.vector))
}

if(conf == TRUE){
#Bootstrap
set.seed(42)
myboot <- boot(matrix.group, function(matrix.group, ind, mod, preds) {
  matrix.group$Value <- fitted(mod) + residuals(mod)[ind]
  tryCatch(predict(nls(Value ~ SSasymp(Step, Asym, resp0, lrc), data = matrix.group), newdata = preds), error = function(e) preds$Step * NA)}, mod = mod, preds = preds, R = 1e4)

#Calculate confidence intervals
CI <- t(sapply(seq_len(nrow(preds)), function(i) boot.ci(myboot, type = "bca", index = i)$bca[4:5]))
colnames(CI) <- c("Mod.lwr","Mod.upr")
group.result <- cbind(raref.obs.NA,Modelled, unique(CI))
}else{
group.result <- cbind(raref.obs.NA[,c(1:2)],Modelled,coverage.vector,completeness.vector)
colnames(group.result)[4] <- c("Coverage")
colnames(group.result)[5] <- c("Completeness")

}

all.results[[g]][["Steps"]] <- group.result
all.results[[g]][["Model.fit"]] <- correlation
all.results[[g]][["Asymptote"]] <- unname(asymp.est)
all.results[[g]][["Coverage"]] <- unname(coverage)
all.results[[g]][["Min.sample"]] <- est.sample.size
all.results[[g]][["Completeness"]] <- completeness

}#close group  

return(all.results)
}

}
