Alpha diversity
| Arguments | Description |
| ------------- | ------------- |
| otutable | A count table (matrix/data.frame) indicating the absolute or relative OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs. |
| qvalue | A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match. |
| weight | A vector indicating the relative weight of each sample. The order needs to be identical to the order of the samples in the OTU table. The values need to sum up to 1. If empty, all samples are weighed the same.  |

````R
alpha.div(otutable=otu.table,qvalue=1,weight=weight.vector)
alpha.div(otutable=otu.table,qvalue=1,tree=tree)
alpha.div(otu.table,1,tree,weight.vector)
````
