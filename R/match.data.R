#' Match data
#' @title Match data
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV tree names
#' @description Filter count tables and OTU/ASV phylogenetic trees to match OTUs/ASVs present in both data files..
#' @param otutable A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @param tree An ultrametic tree of class 'phylo'.
#' @param output Whether to output a filtered OTU table (matrix) or a filtered OTU tree (phylo).
#' @seealso \code{\link{hill.div}}, \code{\link{index.div}}
#' @examples
#' match.data(bat.diet.otutable,bat.diet.tree,output="otutable")
#' match.data(bat.diet.otutable,bat.diet.tree,output="tree")
#' @export

match.data <- function(otutable,tree,output){

#Obtain data
  #OTUs in OTU table
  otutable.otus <- rownames(otutable)
  #Samples in OTU table
  otutable.samples <- colnames(otutable)
  #OTUs in tree
  tree.otus <- tree$tip.label

if(missing(output)){
  if((length(setdiff(tree.otus,otutable.otus)) == 0) & (length(setdiff(otutable.otus,tree.otus)) == 0)){message("OTUs in the OTU table and OTU tree match perfectly.")}
  if((length(setdiff(tree.otus,otutable.otus)) > 0) & (length(setdiff(otutable.otus,tree.otus)) == 0)){message("The OTU tree contains OTUs absent in the OTU table. Filter the OTU tree")}
  if((length(setdiff(tree.otus,otutable.otus)) == 0) & (length(setdiff(otutable.otus,tree.otus)) > 0)){message("The OTU table contains OTUs absent in the OTU tree. Filter the OTU table")}
  if((length(setdiff(tree.otus,otutable.otus)) > 0) & (length(setdiff(otutable.otus,tree.otus)) > 0)){message("The OTU table contains OTUs absent in the OTU tree and the OTU tree contains OTUs absent in the OTU table. Filter both files")}
  output="NA"
}

#Output OTU table
if(output == "otutable"){
  if(length(setdiff(tree.otus,otutable.otus)) > 0){
    message("All OTUs present in the OTU table are present in the tree, but the tree contains OTUs absent in the OTU table. Remember to filter the tree.")
  }
  if(length(setdiff(otutable.otus,tree.otus)) > 0){
    OTUs.to.drop <- setdiff(otutable.otus,tree.otus)
    otutable.filt <- otutable[!(row.names(otutable) %in% OTUs.to.drop), ]
    OTUs.to.drop.string <- paste(OTUs.to.drop,collapse=", ")
    message("The following OTUs were removed from the OTU table for being absent in the tree: ",OTUs.to.drop.string)
    return(otutable.filt)
  }
  if((length(setdiff(tree.otus,otutable.otus)) == 0) & (length(setdiff(otutable.otus,tree.otus)) == 0)){
    message("OTUs in the OTU table and OTU tree match perfectly. No new OTU table was created.")
  }

}

#Output tree
if(output == "tree"){
  if(length(setdiff(tree.otus,otutable.otus)) > 0){
    OTUs.to.drop <- setdiff(tree.otus,otutable.otus)
    tree.filt <- drop.tip(tree,OTUs.to.drop)
    OTUs.to.drop.string <- paste(OTUs.to.drop,collapse=", ")
    message("The following OTUs were removed from the tree for being absent in the OTU table: ",OTUs.to.drop.string)
    return(tree.filt)
  }
  if(length(setdiff(otutable.otus,tree.otus)) > 0){
    message("All OTUs present in the tree are present in the OTU table, but the OTU table contains OTUs absent in the tree. Remember to filter the OTU table.")
  }
  if((length(setdiff(tree.otus,otutable.otus)) == 0) & (length(setdiff(otutable.otus,tree.otus)) == 0)){
    message("OTUs in the OTU tree and OTU table match perfectly. No new tree was created.")
  }
}

}
