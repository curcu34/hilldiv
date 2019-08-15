# div.profile()

Create diversity profiles of a single or multiple samples displayed independently or aggregated in groups. Diversity profiles show the relation between the order of diversity (q-value) and the respective Hill numbers, thus providing information about the richness and evenness of a sample at a glance.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| abund | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |
| qvalues | A vector of sequential orders of diversity (default from 0 to 5). order=seq(from = 0, to = 5, by = (0.1))|
| tree | A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).  |
| values | If 'TRUE' the Hill numbers for each order of diversity for each sample or group (if hierarchy table is provided) are printed on screen. values="FALSE".  |
| hierarchy | A two-column matrix indicating the relation between samples (first column) and groups (second column).  |
| level | Only if hierachy is used. If 'alpha', alpha diversity of the groups is plotted instead of gamma diversity. level="gamma"  |
| colour | The number of vector items (colours, e.g. '#34k235'), must equal the number of samples or groups that are intended to plot.  |
| log | If 'TRUE' the Hill numbers are transformed to the logarithmic scale. This is useful when there are large differences between q values (e.g. sharp drop from q=0 to q=1), which might complicate visualization. log="FALSE"  |

## Examples
````R
# One sample
div.profile(otu.vector)
# Multiple individual samples (first 5 samples of the OTU table)
div.profile(otu.table[,c(1:5)])
# Multiple groups (aggregated samples)
div.profile(otu.table,hierarchy=hierarchy.table,colour=c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7"))````
````

<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.one.png" width="350" title="One sample">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.multiple.png" width="350" title="Multiple samples">
