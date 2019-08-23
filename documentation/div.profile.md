# div.profile()

Create diversity profile vector (single sample or system) or tables (multiple samples or groups) from count tables.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| count | A vector or a matrix indicating the (relative) OTU/ASV counts of one or multiple samples. If a matrix is provided, columns must refer to samples and rows to OTUs. |
| qvalues | A vector of sequential orders of diversity (default from 0 to 5). order=seq(from = 0, to = 5, by = (0.1))|
| tree | A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). Use the function match.data() if the OTU names do not match. |
| values | If 'TRUE' the Hill numbers for each order of diversity for each sample or group (if hierarchy table is provided) are printed on screen. values="FALSE".  |
| hierarchy | A two-column matrix indicating the relation between samples (first column) and groups (second column).  |
| level | Whether to compute alpha or gamma diversities of the system or the groups specified in the hierarchy table.  |

## Examples
````R
data(bat.diet.otutable)
#Diversity profile of a single sample (default profile)
div.profile(bat.diet.otutable[,1])
#Diversity profile of a single sample (custom profile from 0 to 3)
div.profile(bat.diet.otutable[,1],qvalues=c(seq(from = 0, to = 3, by = (0.1)))

#Diversity profile of multiple samples
div.profile(bat.diet.otutable)
#Diversity profile of alpha diversity of multiple samples (entire count table).
div.profile(bat.diet.otutable,level="alpha")
#Diversity profile of gamma diversity of multiple samples (entire count table)
div.profile(bat.diet.otutable,level="gamma")

#Diversity profiles of groups determined by hierarchy tables
data(bat.diet.hierarchy)
#Alpha diversity
div.profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="alpha")
#Gamma diversity
div.profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="gamma")
````
