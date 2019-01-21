div.test <- function(otutable,qvalue,hierarchy,tree){ 
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(hierarchy)) stop("Hierarchy table is necessary to contrast groups of samples")

if(missing(tree)){
div.values <- hill.div(otutable,qvalue)
}else{
if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")  
if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
div.values <- hill.div(otutable,qvalue,tree)
}
div.values.groups <- merge(t(t(div.values)),hierarchy,by.x="row.names",by.y="Sample")
colnames(div.values.groups) <- c("Sample","Value","Group")
  
#Data distribution (normality and homogeneity) assessment
shapiro <- shapiro.test(div.values.groups$Value)
barlett <- bartlett.test(Value ~ Group, data= div.values.groups)
if((shapiro$p.value >= 0.05) & (barlett$p.value >= 0.05)){
norm.homo=TRUE
}else{
norm.homo=FALSE
}
  
#Statistical test
if(length(unique(div.values.groups$Group)) == 2){
    if(norm.homo=TRUE){
    method <- "Student's t-Test"
    test <- t.test(Value ~ Group, data = div.values.groups)  
    results <- list(normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "Student's t-Test",result=c(t=unname(test$statistic),df=unname(test$parameter),p.value=unname(test$p.value)))
    }else{
    test <- wilcox.test(Value ~ Group, data = div.values.groups)
    results <- list(normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "Wilcoxon Rank Sum Test",result=c(W=unname(test$statistic),df=unname(test$parameter),p.value=unname(test$p.value)))
    }
}else{
    if(norm.homo=TRUE){
    test <- summary(aov(Value ~ Group, data = div.values.groups))
    results <- list(normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "ANOVA",result=test)
    }else{
    test <- kruskal.test(Value ~ Group, data = div.values.groups)
    results <- list(normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "Kruskal-Wallis Test",result=c(chi.squared=unname(test$statistic),df=unname(test$parameter),p.value=unname(test$p.value)))
    }
}
 
  
  result = c(t = res[1]), p.value = res[2], m, data.name = DNAME)
  
  groups <- 
  
return(test)
}
