# tree.depth()

Computes phylogenetic tree depth based from a phylogenetic tree and a vector of (relative) abundances.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match. |
| abund | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |

## Examples
````R
tree.depth(tree=otu.tree,abund=otu.table)
tree.depth(otu.tree,otu.table)
````
