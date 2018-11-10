DiverHill is a R package for the analysis of diversity based on Hill numbers

# Functions
### true.div()
True (phylo)diversity computation of individual samples from vectors (one sample) or OTU tables (multiple samples). Using the tree argument yields phylodiversity values. Note that if using a tree the tip labels and the 'names' (vectors) or 'rownames' (matrices) need to be identical.

|Argument|Feature|Type|Notes|
| ------------- | ------------- | ------------- |------------- |
| **abund**  | Mandatory  |vector,matrix,data.frame| A vector or a matrix |
| **q.value**  | Mandatory |>0 |Possitive number, can contain decimals   |
| **tree**   | Optional  |phylo   |Phylogenetic tree of class phylo   |
    
````R
#EXAMPLES#
true.div(vector,0)
true.div(otu.table,1)
true.div(otu.table,qvalue=1)
true.div(otu.table,1,tree)
true.div(otu.table,qvalue=1,tree=tree)
true.div(otu.table[,1],qvalue=1)
9.145646
true.div(otu.table[,c(1:3)],qvalue=1)
    TUL1     TUL2     TUL3 
9.145646 8.686439 7.884177 
````
### div.comp.test()
Diversity comparison between two or multiple groups of samples. If the tree argument is used the test compares phylodiversity values. 
````R
#EXAMPLES#
div.contrast(otu.table,1,hierarchy.table)
div.contrast(otu.table,1,hierarchy.table,tree)
````
### div.comp.plot()
Visual comparison between the diversity levels of two or multiple groups of samples. The chart argument enables selecting between boxplot and jitter plot.
````R
#EXAMPLES#
div.comp.plot(otu.table,1,hierarchy.table)
div.comp.plot(otu.table,1,hierarchy.table,chart="boxplot") #default option
div.comp.plot(otu.table,2,hierarchy.table,chart="jitter") #plot jitter plot instead of boxplot
````

### alpha.div() / alpha.phylodiv()
Alpha (phylo)diversity of a dataset (OTU table) with multiple samples.

### gamma.div() / gamma.phylodiv()
Gamma (phylo)diversity of a dataset (OTU table) with multiple samples.

### div.part() / phylodiv.part()
(Phylo)Diversity partitioning of a dataset (OTU table) with multiple samples and hierarchical levels (2, 3 or 4 levels; e.g. sample > population > species > overall).
````R
#EXAMPLES#
div.part(otu.table,0)
div.part(otu.table,0,hierarchy.table)
div.part(otu.table,qvalue=0,hierarchy=hierarchy.table)
````

###  pair.dis()
Pairwise dissimilarity measurement yielding true beta diversity, homogeneity, overlap and turnover values.
````R
#EXAMPLES#
pair.dis(otu.table,0)
pair.dis(otu.table,qvalue=1)
pair.dis(otu.table,1,measure="overlap")
pair.dis(otu.table,1,measure=c("homogeneity","overlap"))
pair.dis(otu.table,1,measure=c("homogeneity","overlap"),hierarchy.table[,c(1:2)])
````

### div.profile() / phylodiv.profile() - chart
(Phylo)Diversity profiles of individual samples or groups of samples.
````R
#EXAMPLES#
div.profile(otu.table[,1])
div.profile(otu.table[,c(1:5)])
````

<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.one.png" width="350" title="One sample">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.multiple.png" width="350" title="Multiple samples">

# Input data formats

## OTU table
For functions true.div, true.phylodiv, alpha.div, alpha.phylodiv, gamma.div, gamma.phylodiv, div.part, phylodiv.part...

|       | Sample1 | Sample2 |Sample3  |SampleN  |
| ------------- | ------------- | ------------- | ------------- |------------- |
| OTU1  | 0.604   |0.543    |0.204    |...    |
| OTU2  | 0.050   |0.210    |0.450    |...  |
| OTU3  | 0.153   |0.105    |0.314    |...  |
| OTUN  | ...   |...    |...    |...  |

The sum of relative abundances of all samples must be 1. Check that by using colSums(otutable)

## Hierarchy table
For functions div.part and phylodiv.part

| Sample | Group |Subgroup  |
| ------------- | ------------- |------------- |
| Sample1  | Group1   |Subgroup1  |
| Sample2  | Group1   |Subgroup2  |
| Sample3  | Group1   |Subgroup2  |
| Sample4  | Group2   |Subgroup3  |
| Sample5  | Group3   |Subgroup4  |
| ...  | ...   |...  |
| SampleN  | GroupN   |SubgroupN  |

* Sample names and OTU table row names must match
* Groups and subgroups must be nested
