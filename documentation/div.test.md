# div.test()

Diversity comparison test between groups of samples.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| otutable | A beta diversity value based on Hill numbers. |
| qvalue |  A positive integer or decimal number (>=0), usually between 0 and 3. |
| hierarchy | A two-column matrix indicating the relation between samples (first column) and groups (second column). |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function [match.data()](documentation/match.data.md) if the OTU names do not match.  |

## Examples
````R
div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
div.test(otu.table,qvalue=1,hierarchy=hierarchy.table,tree=tree)
div.test(otu.table,2,hierarchy.table,tree)
````
