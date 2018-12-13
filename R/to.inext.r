to.inext <- function(otutable,hierarchy,type){

if(missing(type)){type="incidence_raw"}
	
#Without hierarchy (single list)
if(missing(hierarchy)){
if(type == "incidence_raw"){
otutable[which(otutable != 0)] <- 1
result <- otutable}
if(type == "incidence_freq"){
result <- c(ncol(otutable),rev(sort(rowSums(otutable != 0))))}
result <- result[result > 0]
return(result)
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
	if(type == "incidence_raw"){
	otutable.subset[which(otutable.subset != 0)] <- 1}
	if(type == "incidence_freq"){
	list.subset <- c(ncol(otutable.subset),rev(sort(rowSums(otutable.subset != 0))))}
	results[[g]] <- list.subset
}
}
return(results)
    
}
