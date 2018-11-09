DiverHill is a R package for the analysis of diversity based on Hill numbers

## Functions
#### true.div() / true.phylodiv()
True (phylo)diversity computation of individual samples from vectors (1 sample) or OTU tables (multiple samples).
````R
#EXAMPLES#
true.div(vector,0)
true.div(otu.table,1)
true.div(otu.table,qvalue=1)
true.phylodiv(otu.table,1,tree)
true.phylodiv(otu.table,qvalue=1,tree=tree)
true.div(otu.table[,1],qvalue=1)
9.145646
true.div(otu.table[,c(1:3)],qvalue=1)
    TUL1     TUL2     TUL3 
9.145646 8.686439 7.884177 
````

#### alpha.div() / alpha.phylodiv()
Alpha (phylo)diversity of a dataset (OTU table) with multiple samples.

#### gamma.div() / gamma.phylodiv()
Gamma (phylo)diversity of a dataset (OTU table) with multiple samples.

#### div.part() / phylodiv.part()
(Phylo)Diversity partitioning of a dataset (OTU table) with multiple samples and hierarchical levels (2, 3 or 4 levels; e.g. sample > population > species > overall).

####  pair.dis()
Pairwise dissimilarity measurement yielding true beta diversity, homogeneity, overlap and turnover values.
````R
#EXAMPLES#
pair.dis(otu.table,0)
pair.dis(otu.table,qvalue=1)
pair.dis(otu.table,1,measure="overlap")
pair.dis(otu.table,1,measure=c("homogeneity","overlap"))
````

#### div.profile() / phylodiv.profile() - chart
(Phylo)Diversity profiles of individual samples or groups of samples.

## Input data formats

### OTU table
For functions true.div, true.phylodiv, alpha.div, alpha.phylodiv, gamma.div, gamma.phylodiv, div.part, phylodiv.part...

|       | Sample1 | Sample2 |Sample3  |SampleN  |
| ------------- | ------------- | ------------- | ------------- |------------- |
| OTU1  | 0.604   |0.543    |0.204    |...    |
| OTU2  | 0.050   |0.210    |0.450    |...  |
| OTU3  | 0.153   |0.105    |0.314    |...  |
| OTUN  | ...   |...    |...    |...  |

The sum of relative abundances of all samples must be 1. Check that by using colSums(otutable)

### Hierarchy table
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
