to.inext <- function(otutable,hierarchy){

#Without hierarchy (single list)
if(missing(hierarchy)){
list <- c(ncol(otutable.subset),rev(sort(rowSums(otutable.subset != 0))))
return(results)
}

#With hierarchy (multiple lists)
if(!missing(hierarchy)){
colnames(hierarchy) <- c("Sample","Group")
groups <- sort(unique(hierarchy$Group))   

results <- list()
for (g in groups){
	samples <- as.character(hierarchy[which(hierarchy$Group == g),1])
	otutable.subset <- otutable[,samples]
	otutable.subset <- as.data.frame(otutable.subset[apply(otutable.subset, 1, function(z) !all(z==0)),])
	list.subset <- c(ncol(otutable.subset),rev(sort(rowSums(otutable.subset != 0))))
	results[[g]] <- list.subset
}
}
return(results)
    
}
