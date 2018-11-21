gamma.curve <- function(otutable,qvalue,tree,hierarchy,steps,iter,summary){
	
if(missing(steps)){steps=100}
if(missing(iter)){iter=50}
if(missing(summary)){summary=FALSE}

#################################
# IF HIERARCHY IS NOT SPECIFIED #
#################################

if(missing(hierarchy) == TRUE){

#Define step vector
max.steps <- ncol(otutable)
selected.steps <- max.steps * steps / 100
interval <- max.steps/selected.steps
step.vector <- seq(from = 1, to = max.steps, by = interval)
if(step.vector[length(step.vector)] != max.steps){ step.vector <- c(step.vector,max.steps)} #Add last number if is not included	

#Loop across steps and number of iterations	
matrix <- c()
for (i in 1:iter){#open iteration
  vector <- c()
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

}#close no-hierarchy function

#############################
# IF HIERARCHY IS SPECIFIED #
#############################		
		
if(ncol(hierarchy) == 2){

#Define step vector	
max.steps <- max(table(hierarchy[,2]))
selected.steps <- max.steps * steps / 100
interval <- max.steps/selected.steps
step.vector <- seq(from = 1, to = max.steps, by = interval)
if(step.vector[length(step.vector)] != max.steps){ step.vector <- c(step.vector,max.steps)} #Add last number if is not included	

#Define groups	
groups <- as.character(unique(hierarchy[,2]))
	
#Loop across steps, groups and number of iterations	
matrix <- c()
for (g in groups){#open group
samples <- hierarchy[which(hierarchy[,2] == g),1]
otutable.group <- otutable[,samples]
    for (i in 1:iter){#open iteration
	  vector <- c()
	  for (s in step.vector){#open steps
	      subset <- sample(max.steps, s, replace = FALSE)
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
	    vector <- c(vector,value)
	    }#close steps
	 vector <- c(g,vector)
	 matrix <- rbind(matrix,vector)    
    }#close iteration
}#close group  

#Generate summary table	
rownames(matrix) <- c(1:nrow(matrix))	
colnames(matrix) <- c("group",step.vector)	
matrix <- as.data.frame(matrix)	
matrix[,-1] = apply(matrix[,-1], 2, function(x) as.numeric(as.character(x)));
mean <- aggregate(matrix[,-1],by=list(matrix[,1]),FUN=mean)
mean.melt <- melt(mean)
stderr <- aggregate(matrix[,-1],by=list(matrix[,1]),function(x) sd(x)/sqrt(length(x)))
stderr.melt <- melt(stderr)
summary.table <- cbind(mean.melt,stderr.melt[,3])
colnames(summary.table) <- c("Group","step","mean","stderr")		    

#Generate curve plot
getPalette = colorRampPalette(brewer.pal(length(groups), "Paired"))
curve.plot <- ggplot(summary.table, aes(x = step, y = mean, group=Group)) +
       geom_line(aes(colour=Group)) +
       geom_ribbon(aes(ymin = mean - stderr, ymax = mean + stderr, fill=Group), alpha = 0.2) +
       xlab("Sample size") + 
       ylab(if(missing(tree)){"Effective number of OTUs"}else{"Effective number of lineages"}) +
       scale_colour_manual(values = getPalette(length(groups))) + 
       scale_fill_manual(values = getPalette(length(groups))) + 
       theme_minimal()
print(curve.plot)

if(summary == TRUE){
return(summary.table)
}
		    
}#close hierarchy function
                
}#close function
