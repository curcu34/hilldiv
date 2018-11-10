div.comp.test <- function(otutable,qvalue,hierarchy){ 
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(hierarchy)) stop("Hierarchy table is necessary to contrast groups of samples")

div.values <- true.div(otutable,qvalue)
div.values.groups <- merge(t(t(div.values)),hierarchy,by.x="row.names",by.y="Sample")
colnames(div.values.groups) <- c("Sample","Value","Group")

#Statistical test
if(length(unique(div.values.groups$Group)) == 2){
test <- wilcox.test(Value ~ Group, data = div.values.groups)
}else{
test <- kruskal.test(Value ~ Group, data = div.values.groups)
}
return(test)
}
