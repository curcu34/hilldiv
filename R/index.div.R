index.div <- function(abund,tree,index){

#index.div(bat.diet.otutable)
#index.div(bat.diet.otutable,index="shannon")
#index.div(bat.diet.otutable,tree=bat.diet.tree,index="faith")

#Data input control
if(missing(abund)) stop("Abundance data is missing")
if(missing(index) & missing(tree)){index="richness"}
if(missing(index) & !missing(tree)){index="faith"}

#Input data type identification
if(is.null(dim(abund)) == TRUE){
inputtype="onesample"
}else{
inputtype="multiplesamples"
}

#INDICES

#Richness
if(index == "richness"){
  if(!missing(tree)) warning("Phylogenetic tree was not used to compute richness.")
  #Function
  richness <- function(vector){sum(vector != 0)}

  #One sample
  if(inputtype == "onesample"){
  div <- richness(abund)
  return(div)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  divs <- apply(abund, 2, function(x) richness(x))
  return(divs)
  }
}

#Shannon index
if(index == "shannon"){
  if(!missing(tree)) warning("Phylogenetic tree was not used to compute Shannon index.")
  #Function
  shannon <- function(vector){
    pi <- tss(vector[vector != 0])
    div <- -sum(pi*log(pi))
    }

  #One sample
  if(inputtype == "onesample"){
  div <- shannon(abund)
  return(div)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  divs <- apply(abund, 2, function(x) shannon(x))
  return(divs)
  }

}

#Simpson index
if(index == "simpson"){
  if(!missing(tree)) warning("Phylogenetic tree was not used to compute Simpson index.")
  #Function
  simpson <- function(vector){
    pi <- tss(vector[vector != 0])
    div <- 1 - (sum(pi^2))
    }

  #One sample
  if(inputtype == "onesample"){
  div <- simpson(abund)
  return(div)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  divs <- apply(abund, 2, function(x) simpson(x))
  return(divs)
  }
}

#Faith's PD
if(index == "faith"){
if(missing(tree)) stop("Faith's PD cannot be computed without a phylogenetic tree")
  #Function
  faith <- function(vector,tree){
  Li <- tree$edge.length
  ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
  ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
  Li <- Li[ai != 0]
  phylodiv <- sum(Li)
  }

  #One sample
  if(inputtype == "onesample"){
  div <- faith(abund,tree)
  return(div)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  divs <- apply(abund, 2, function(x) faith(x,tree))
  return(divs)
  }

}

#Allen's H
if(index == "allen"){
if(missing(tree)) stop("Allen's H cannot be computed without a phylogenetic tree")
  #Function
  allen <- function(vector,tree){
  vector <- tss(vector)
  Li <- tree$edge.length
  ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
  ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
  T <- sum(Li * ai)
  Li <- Li[ai != 0]
  ai <- ai[ai != 0]
  phylodiv <- -sum(Li*ai*log(ai))
  }

  #One sample
  if(inputtype == "onesample"){
  phylodiv <- allen(abund,tree)
  return(phylodiv)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  phylodivs <- apply(abund, 2, function(x) allen(x,tree))
  return(phylodivs)
  }

}

#Rao's Q
if(index == "rao"){
if(missing(tree)) stop("Rao's Q cannot be computed without a phylogenetic tree")
  #Function
  rao <- function(vector,tree){
  vector <- tss(vector)
  pij <- outer(vector, vector)
  phylodist <- cophenetic.phylo(tree)/2
  phylodist <- phylodist[rownames(pij),colnames(pij)]
  phylodiv <- sum(phylodist * pij)
  }

  #One sample
  if(inputtype == "onesample"){
  phylodiv <- rao(abund,tree)
  return(phylodiv)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  phylodivs <- apply(abund, 2, function(x) rao(x,tree))
  return(phylodivs)
  }

}

}
