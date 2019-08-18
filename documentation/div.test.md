# div.test()

Diversity comparison test between groups of samples. The function automatically assesses whether the data meets the properties for parametric statistics and performs the appropriate test accordingly: Students' T, ANOVA, Wilcoxon or Kruskal-Wallis. If the posthoc argument is set as TRUE, multiple group comparisons are complemented with post hoc pairwise tests, either Tukey test (parametric) or Dunn test with Benjamini-Hochberg correction (non-parametric). 

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| otutable | A beta diversity value based on Hill numbers. |
| qvalue |  A positive integer or decimal number (>=0), usually between 0 and 3. |
| hierarchy | A two-column matrix indicating the relation between samples (first column) and groups (second column). |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function [match.data()](match.data.md) if the OTU names do not match.  |
| posthoc | Whether to run post hoc pairwise analyses or not. If TRUE, an ANOVA will be complemented with a Tukey test and a Kruskal-Wallis test will be complemented with a Dunn test. |

## Examples
````R
div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
div.test(otu.table,qvalue=1,hierarchy=hierarchy.table,tree=tree)
div.test(otu.table,2,hierarchy.table,tree)
div.test(otu.table,qvalue=1,hierarchy=hierarchy.table,posthoc=TRUE)
````
