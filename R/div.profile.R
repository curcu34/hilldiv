#' Diversity profile
#' @title Diversity profile
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords alpha gamma beta hill
#' @description Create diversity profiles of a single or multiple samples displayed independently or aggregated in groups.
#' @param abund A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @param qvalues A vector of sequential orders of diversity (default from 0 to 5). order=seq(from = 0, to = 5, by = (0.1))
#' @param tree A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).
#' @param values If 'TRUE' the Hill numbers for each order of diversity for each sample or group (if hierarchy table is provided) are printed on screen. values="FALSE".
#' @param hierarchy A two-column matrix indicating the relation between samples (first column) and groups (second column).
#' @param level Only if hierachy is used. If 'alpha', alpha diversity of the groups is plotted instead of gamma diversity. level="gamma"
#' @param colour The number of vector items (colours, e.g. '#34k235'), must equal the number of samples or groups that are intended to plot.
#' @param log If 'TRUE' the Hill numbers are transformed to the logarithmic scale. This is useful when there are large differences between q values (e.g. sharp drop from q=0 to q=1), which might complicate visualization. log="FALSE"
#' @return A figure or a table of sequential diversity values.
#' @seealso \code{\link{hill.div}}, \code{\link{div.part}}
#' @examples
#' div.profile(abund=otu.table[,1],qvalues=seq(from = 0, to = 5, by = (0.1)),tree=tree,values="TRUE")
#' div.profile(vector)
#' div.profile(otu.table,hierarchy=hierarchy.table)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources. Early view.\cr\cr
#' Chao, A., Chiu, C.‐H., & Jost, L. (2014). Unifying species diversity, phylo‐ genetic diversity, functional diversity, and related similarity and dif‐ ferentiation measures through hill numbers. Annual Review of Ecology Evolution and Systematics, 45, 297–324.
#' @export

div.profile <- function(abund,qvalues,tree,values,hierarchy,level,colour,log){

#Quality-check and warnings
if(missing(abund)) stop("The abundance data is missing")
if(missing(qvalues)) {qvalues= seq(from = 0, to = 5, by = (0.1))}
if(missing(values)) {values= "FALSE"}
if(missing(level)) {level= "gamma"}
if(missing(log)) {log= "FALSE"}

#If input data is a vector
if(is.null(dim(abund)) == TRUE){
    profile <- c()
    for (o in qvalues){
    if(missing(tree)){
        div.value <- hilldiv::hill.div(abund,o)
        }else{
        div.value <- hilldiv::hill.div(abund,o,tree)
    }
    profile <- c(profile,div.value)
    }
    profile.melted <- as.data.frame(cbind(qvalues,profile))
    colnames(profile.melted) <- c("Order","Profile")
    plot <- ggplot(profile.melted , aes(x = Order, y = Profile)) +
           geom_line() +
           xlab("Order of diversity") + ylab("Effective number of OTUs") +
           theme_minimal()
    print(plot)
    if(values == "TRUE"){
    return(profile)
    }
}

#If input data is an OTU table
if(is.null(dim(abund)) == FALSE){

    if(dim(abund)[1] < 2) stop("The OTU table only less than 2 OTUs")
    if(dim(abund)[2] < 2) stop("The OTU table contains less than 2 samples")

    profile <- c()
    if(missing(hierarchy)){
        for (o in qvalues){
            if(missing(tree)){
            div.values <- hilldiv::hill.div(abund,o)
            }else{
            div.values <- hilldiv::hill.div(abund,o,tree)
            }
        profile <- rbind(profile,div.values)
        }
        rownames(profile) <- qvalues
        profile.melted <- as.data.frame(melt(profile))
        colnames(profile.melted) <- c("Order","Sample","Value")
        profile.melted[,1] <- as.numeric(as.character(profile.melted[,1]))
        profile.melted[,3] <- as.numeric(as.character(profile.melted[,3]))
        if(log == "TRUE"){profile.melted[,3] <- log(profile.melted[,3])}

        #Declare colours
	if(missing(colour) || (length(colour) != ncol(abund))){
        getPalette <- colorRampPalette(brewer.pal(ncol(abund), "Paired"))
        colour <- getPalette(ncol(abund))
        }

        #Plot
        plot <- ggplot(profile.melted , aes(x = Order, y = Value, group=Sample, colour=Sample)) +
        geom_line() +
        xlab("Order of diversity") +
        ylab(if(log == "TRUE"){"Effective number of OTUs (log-transformed)" }else{"Effective number of OTUs"}) +
        scale_colour_manual(values = colour) +
        theme_minimal()
        print(plot)

    }else{
    colnames(hierarchy) <- c("Sample","Group")
    groups <- sort(unique(hierarchy$Group))
        for (g in groups){
            samples <- as.character(hierarchy[which(hierarchy$Group == g),1])
            abund.subset <- abund[,samples]
            abund.subset <- as.data.frame(abund.subset[apply(abund.subset, 1, function(z) !all(z==0)),])
            if(!missing(tree)){
            missing.otus <- setdiff(tree$tip.label,rownames(abund.subset))
            tree.subset <- drop.tip(tree,missing.otus)
            }
                 for (o in qvalues){
                        if(missing(tree)){
                            if(level == "gamma"){div.value <- hilldiv::gamma.div(abund.subset,o)}
                            if(level == "alpha"){div.value <- hilldiv::alpha.div(abund.subset,o)}
                            if(level == "incidence"){div.value <- hilldiv::hill.div(rowSums(abund.subset != 0)/sum(rowSums(abund.subset != 0)),o)}
                        }else{
                            if(level == "gamma"){div.value <- hilldiv::gamma.div(abund.subset,o,tree.subset)}
                            if(level == "alpha"){div.value <- hilldiv::alpha.div(abund.subset,o,tree.subset)}
                            if(level == "incidence"){div.value <- hilldiv::hill.div(rowSums(abund.subset != 0)/sum(rowSums(abund.subset != 0)),o,tree.subset)}
                        }
                 profile <- rbind(profile,cbind(g,div.value))
                 }
          }

    profile <- as.data.frame(cbind(profile,rep(qvalues,length(groups))))
    profile[,2] <- as.numeric(as.character(profile[,2]))
    profile[,3] <- as.numeric(as.character(profile[,3]))
    colnames(profile) <- c("Group","Value","Order")
    if(log == "TRUE"){profile[,2] <- log(profile[,2])}

    #Declare colours
    if(missing(colour) || (length(colour) != length(groups))){
    getPalette <- colorRampPalette(brewer.pal(length(groups), "Paired"))
    colour <- getPalette(length(groups))
    }

    #Plot
    plot <- ggplot(profile , aes(x = Order, y = Value, group=Group, colour=Group)) +
    geom_line() +
    xlab("Order of diversity") +
    ylab(if((log == "TRUE") & missing(tree)){"Effective number of OTUs (log-transformed)"}else if((log == "TRUE") & !missing(tree)){"Effective number of lineages (log-transformed)"}else if((log == "FALSE") & !missing(tree)){"Effective number of lineages"}else{"Effective number of OTUs"}) +
    scale_colour_manual(values = colour) +
    theme_minimal()
    print(plot)

    }




    if(values == "TRUE"){
    return(profile)
    }


}

}
