# div.profile.plot()

Create diversity profiles of a single or multiple samples displayed independently or aggregated in groups. Diversity profiles show the relation between the order of diversity (q-value) and the respective Hill numbers, thus providing information about the richness and evenness of a sample at a glance.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| profile | A div.profile() object or a vector/matrix containg diversity profile(s), with columns indicating samples/groups and rows indicating orders of diversity (q-values). |
| colour | A vector of RGB colours e.g. c("#34k235","#99cc00"). The number of vector items, must equal the number of samples or groups that are intended to plot. |
| log | Whether to transform Hill numbers to logarithmic scale (TRUE) or not (FALSE). This is useful when there are large differences between q values (e.g. sharp drop from q=0 to q=1), which might complicate visualization. Default: log="FALSE"  |
| legend | Whether to display the legend (TRUE) or not (FALSE) in diversity profiles containing multiple samples/groups. Default TRUE in multi-sample charts.  |

## Examples
````R
# Multiple samples
profile <- div.profile(bat.diet.otutable)
div.profile.plot(profile)
# One system
profile <- div.profile(bat.diet.otutable,level="gamma")
div.profile.plot(profile)
# Multiple groups
profile <- div.profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy)
div.profile.plot(profile,colour=c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7"),log=TRUE)
````

<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.one.png" width="350" title="One sample">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.multiple.png" width="350" title="Multiple samples">
