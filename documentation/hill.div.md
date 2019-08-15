# hill.div

Compute neutral or phylogenetic Hill numbers from a single sample (vector) or count table (matrix). Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs, i.e. the number of equally abundant OTUs that would be needed to give the same value of diversity (Hill 1973; Jost 2006). When all OTUs in a system have the same relative abundances the effective number of OTUs for all q values equals the actual number of OTUs, namely richness. When the relative abundances of the types vary however, then the effective number of OTUs for q>0 values decreases. The higher the heterogeneity between OTUs, the lower the effective number of OTUs. In extreme cases in which the system is dominated by a few equally abundant OTUs, the effective number of OTUs will approach the number of those abundant OTUs. The sensitivity towards abundant and rare OTUs can be modulated using the scaling parameter q, known as the “order” of diversity (Jost 2006). The larger the q value, the higher the importance attributed to abundant OTUs. Three q values are particularly relevant, both for their significance, and their close relationship to popular diversity indices: q=0, q=1 and q=2. When a diversity of order zero (q=0) is applied to the formula, it becomes insensitive to OTU frequencies, thus yielding a richness value. As the relative abundances of OTUs are overlooked, rare OTUs are overweighed. A q value of 1 (in practical terms its limit, e.g. 0.9999, as the Hill number is undefined for q=1) is the value that weighs OTUs by their frequency, without disproportionately favoring either rare or abundant ones (Jost 2006). The value it yields is exactly the exponential of the Shannon index. In fact, q values under unity favor rare OTUs, while values above one favor abundant OTUs. When a q value of 2 is applied, abundant OTUs are overweighed, and the formula yields the multiplicative inverse of the Simpson index. Indeed, common diversity indices can be transformed to Hill numbers (also known as numbers equivalents or true diversities sensu Jost 2006), by applying simple mathematical transformations

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
