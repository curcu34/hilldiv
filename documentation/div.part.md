# div.part()

Diversity partitioning based on Hill numbers. the Hill numbers framework enables the diversity of a system to be partitioned following the multiplicative definition (Jost, 2007; Chao, Chiu, & Hsieh, 2012). Alpha diversity is obtained by computing the Hill numbers from the averaged basic sums of the samples, while gamma diversity is obtained by taking the average of OTU relative abundances across samples, and then computing the Hill numbers of the pooled system. The division of gamma diversity by alpha diversity yields the beta diversity, which quantifies how many times richer an entire system is in effective OTUs (gamma diversity) than its constituent samples are on average (alpha diversity). However, the Hill numbers beta diversity can also be considered an actual diversity value, as the same metric also measures the effective number of equally large and completely distinct samples in a system.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| otutable | An OTU table (matrix/data.frame) indicating the absolute or relative OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs. |
| qvalue | A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| hierarchy | A two-column matrix indicating the relation between samples (first column) and groups (second column).  |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match.  |
| type | Either "abundance" for abundance-based Hill number computation or "incidence" for incidence-based Hill number computation. Incidence-based Hill number computation requires a hierarchy table. |

## Examples
````R
div.part(otu.table,qvalue=1)
div.part(otu.table,qvalue=1,type="abundance")
div.part(otu.table,qvalue=0,tree=tree)
div.part(otu.table,qvalue=0,type="incidence",hierarchy=hierarchy.table)
````
