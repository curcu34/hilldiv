# match.data()

Auxiliary function to match OTUs present in count tables and OTU/ASV phylogenetic trees and ensure phylogenetic Hill numbers are correctly computed.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| otutable | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |
| tree | An ultrametic tree of class 'phylo'. |
| output | Whether to output a filtered OTU table (matrix) or a filtered OTU tree (phylo).  |

## Examples
````R
match.data(bat.diet.otutable,bat.diet.tree,output="otutable")
match.data(bat.diet.otutable,bat.diet.tree,output="tree")
````
