# index.div()
Compute neutral (richness, Shannon index, Simpson index) or phylogenetic (Faith's PD, Allen's H, Rao's Q) diversity indices related to Hill numbers from a vector object (one sample) or an OTU table (matrix or data.frame object; multiple samples). A ultrametric tree object (phylo) is necessary to compute phylogenetic diversity indices. Note that if using a tree the tip labels and the 'names' (vectors) or 'rownames' (matrices) need to be identical. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. If the vector or the OTU table columns do not sum to 1, the data is TSS-normalised.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| abund | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match. |
| index | Diversity index to be computed ("richness", "shannon", "simpson", "faith", "allen", "rao"). Default without tree argument: index="richness". Default with tree argument: index="faith".  |

## Examples
````R
index.div(otu.vector)
index.div(otu.table,index="shannon")
index.div(otu.table,tree=tree)
index.div(otu.table,tree=tree,index="allen")
index.div(otu.table,tree,"rao")
````
