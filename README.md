DiverHill is a R package for the analysis of diversity based on Hill numbers

## Functions
#### true.div() / true.phylodiv()
True (phylo)diversity computation of individual samples from vectors or OTU tables.

#### alpha.div() / alpha.phylodiv()
Alpha (phylo)diversity of a OTU table with multiple samples.

#### gamma.div() / gamma.phylodiv()
Gamma (phylo)diversity of a OTU table with multiple samples.

#### div() / gamma.phylodiv()
(Phylo)Diversity partitioning of a OTU table with multiple samples and hierarchical levels.

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
