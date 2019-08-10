match.data <- function(otutable,tree,output){

#Obtain data
  #OTUs in OTU table
  otutable.otus <- rownames(otutable)
  #Samples in OTU table
  otutable.samples <- colnames(otutable)
  #OTUs in tree
  tree.otus <- tree$tip.label

#Output OTU table
if(output == "otutable"){
  if(length(setdiff(tree.otus,otutable.otus)) > 0){
    warning("All OTUs present in the OTU table are present in the tree, but the tree contains OTUs absent in the OTU table. Remember to filter the tree.")
  }
  if(length(setdiff(otutable.otus,tree.otus)) > 0){
    OTUs.to.drop <- setdiff(otutable.otus,tree.otus)
    otutable.filt <- otutable[!(row.names(otutable) %in% OTUs.to.drop), ]
    OTUs.to.drop.string <- paste(OTUs.to.drop,collapse=", ")
    warning("The following OTUs were removed from the OTU table for being absent in the tree: ",OTUs.to.drop.string)
    return(otutable.filt)
  }
  if((length(setdiff(tree.otus,otutable.otus)) == 0) & (length(setdiff(otutable.otus,tree.otus)) == 0)){
    stop("OTUs in the OTU table and OTU tree match perfectly. No new OTU table was created.")
  }

}

#Output tree
if(output == "tree"){
  if(length(setdiff(tree.otus,otutable.otus)) > 0){
    OTUs.to.drop <- setdiff(tree.otus,otutable.otus)
    tree.filt <- drop.tip(tree,OTUs.to.drop)
    OTUs.to.drop.string <- paste(OTUs.to.drop,collapse=", ")
    warning("The following OTUs were removed from the tree for being absent in the OTU table: ",OTUs.to.drop.string)
    return(tree.filt)
  }
  if(length(setdiff(otutable.otus,tree.otus)) > 0){
    warning("All OTUs present in the tree are present in the OTU table, but the OTU table contains OTUs absent in the tree. Remember to filter the OTU table.")
  }
  if((length(setdiff(tree.otus,otutable.otus)) == 0) & (length(setdiff(otutable.otus,tree.otus)) == 0)){
    stop("OTUs in the OTU tree and OTU table match perfectly. No new tree was created.")
  }
}

}
