gamma.curve <- function(otutable,qvalue,tree,steps,iter,summary){
	
if(missing(steps)){steps=100}
if(missing(iter)){iter=50}
if(missing(summary)){summary=FALSE}
  
matrix <- c()
for (i in 1:iter){#open iteration
  vector <- c()
  max.steps <- ncol(otutable)
  selected.steps <- max.steps * steps / 100
  interval <- max.steps/selected.steps
  step.vector <- seq(from = 1, to = max.steps, by = interval)
  if(step.vector[length(step.vector)] != max.steps){ step.vector <- c(step.vector,max.steps)} #Add last number if is not included
    for (s in step.vector){#open steps
      subset <- sample(max.steps, s, replace = FALSE)
      otutable.subset <- otutable[,subset]
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
    vector <- c(vector,value)
    }#close steps
matrix <- rbind(matrix,vector)
}#close iteration

#Generate summary table
mean <- apply(matrix, 2, mean)
stderr <- apply(matrix, 2, function(x) sd(x)/sqrt(length(x))) 
summary.table <- as.data.frame(cbind(step.vector,mean,stderr))

#Generate curve plot
curve.plot <- ggplot(summary.table, aes(x = step.vector, y = mean)) +
       geom_line() +
       geom_ribbon(aes(ymin = mean - stderr, ymax = mean + stderr), alpha = 0.2) +
       xlab("Sample size") + 
       ylab(if(missing(tree)){"Effective number of OTUs"}else{"Effective number of lineages"}) +
       theme_minimal()

print(curve.plot)

if(summary == TRUE){
return(summary.table)
}
                
}#close function
