# hill.div()

Compute neutral or phylogenetic Hill numbers from a single sample (vector) or count table (matrix). Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs, i.e. the number of equally abundant OTUs that would be needed to give the same value of diversity. 

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| abund | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |
| qvalue | A positive integer or decimal number (>=0), usually between 0 and 3. |
| tree | An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). |
| type | Either "abundance" for abundance-based Hill number computation or "incidence" for incidence-based Hill number computation. In the latter case, the 'abund' object needs to be an OTU table, as incidence-based Hill number computation requires multiple samples. Default: type="abundance".  |

## Examples

hill.div(otu.vector,0)
hill.div(otu.table,1)
hill.div(otu.table,qvalue=2)
hill.div(otu.table,1,tree)
hill.div(otu.table,1,tree,type="incidence")
hill.div(otu.table,qvalue=2,tree=tree)
