**hilldiv** is an R package for the analysis of diversity based on Hill numbers.
# Installation
To install **hilldiv** in your R environment, you need to 1) install devtools, 2) load devtools library, 3) install **hilldiv** using devtools and 4) finally load **hilldiv** library to your environment.

1) install.packages("devtools")
2) library(devtools)
3) install_github("anttonalberdi/hilldiv")
4) library(hilldiv)

# Functions
## Diversity measurement and visualisation
### hill.div()
True (phylo)diversity computation of individual samples from vectors (one sample) or OTU tables (multiple samples). Using the tree argument yields phylodiversity values. Note that if using a tree the tip labels and the 'names' (vectors) or 'rownames' (matrices) need to be identical. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. 

|Argument|Feature|Type|Notes|
| ------------- | ------------- | ------------- |------------- |
| **abund**  | Mandatory  |vector,matrix,data.frame| A vector or a matrix (rows as OTUs)|
| **q.value**  | Mandatory |>0 |Possitive number, can contain decimals   |
| **tree**   | Optional  |phylo   |Phylogenetic tree of class phylo   |
    
````R
#EXAMPLES#
hill.div(vector,0)
hill.div(otu.table,1)
hill.div(otu.table,qvalue=1)
hill.div(otu.table,1,tree)
hill.div(otu.table,qvalue=1,tree=tree)
hill.div(otu.table[,1],qvalue=1)
9.145646
hill.div(otu.table[,c(1:3)],qvalue=1)
    TUL1     TUL2     TUL3 
9.145646 8.686439 7.884177 
````
### div.profile() - chart
(Phylo)Diversity profiles of individual samples or groups of samples.
````R
#EXAMPLES#
div.profile(otu.table[,1])
div.profile(otu.table[,c(1:5)])
````

<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.one.png" width="350" title="One sample">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.multiple.png" width="350" title="Multiple samples">

## Diversity comparison
### div.comp.test()
Diversity comparison between two or multiple groups of samples. If the tree argument is used the test compares phylodiversity values. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. 

|Argument|Feature|Type|Notes|
| ------------- | ------------- | ------------- |------------- |
| **abund**  | Mandatory  |vector,matrix,data.frame| A vector or a matrix (rows as OTUs) |
| **q.value**  | Mandatory |>0 |Possitive number, can contain decimals   |
| **hierarchy.table**  | Mandatory |matrix,data.frame |First column containing sample names, second group names  |
| **tree**   | Optional  |phylo   |Phylogenetic tree of class phylo   |
````R
#EXAMPLES#
div.contrast(otu.table,1,hierarchy.table)
div.contrast(otu.table,1,hierarchy.table,tree)
````
### div.comp.plot() - chart
Visual comparison between the diversity levels of two or multiple groups of samples. The chart argument enables selecting between boxplot and jitter plot.

|Argument|Feature|Type|Notes|
| ------------- | ------------- | ------------- |------------- |
| **abund**  | Mandatory  |vector,matrix,data.frame| A vector or a matrix (rows as OTUs) |
| **q.value**  | Mandatory |>0 |Possitive number, can contain decimals   |
| **hierarchy.table**  | Mandatory |matrix,data.frame |First column containing sample names, second group names  |
| **tree**   | Optional  |phylo   |Phylogenetic tree of class phylo   |
| **chart**   | Optional  |character   |Either "boxplot" or "jitter"   |

````R
#EXAMPLES#
div.comp.plot(otu.table,1,hierarchy.table)
div.comp.plot(otu.table,1,hierarchy.table,chart="boxplot")
div.comp.plot(otu.table,2,hierarchy.table,tree,chart="jitter") 
````
## Diversity estimation from incomplete data sets
### hill.intext()
Interpolation and extrapolation of OTU table diversities based on Hill numbers. Wrapper of iNEXT and iNextPD. 

## Diversity partitioning
### alpha.div()
Alpha (phylo)diversity of a dataset (OTU table) with multiple samples.

### gamma.div()
Gamma (phylo)diversity of a dataset (OTU table) with multiple samples.

### div.part()
(Phylo)Diversity partitioning of a dataset (OTU table) with multiple samples and hierarchical levels (2, 3 or 4 levels; e.g. sample > population > species > overall).
````R
#EXAMPLES#
div.part(otu.table,0)
div.part(otu.table,0,hierarchy.table)
div.part(otu.table,qvalue=0,hierarchy=hierarchy.table)
````

## Dissimilarity measurement and visualisation
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

| Sample | Group |Supergroup  |
| ------------- | ------------- |------------- |
| Sample1  | Group1   |Supergroup1  |
| Sample2  | Group1   |Supergroup1  |
| Sample3  | Group2   |Supergroup2  |
| Sample4  | Group2   |Supergroup2  |
| Sample5  | Group3   |Supergroup2  |
| ...  | ...   |...  |
| SampleN  | GroupN   |SpergroupN  |

* Sample names and OTU table row names must match
* Groups and supergroups must be nested (all samples in a group must be included in the same supergroup)
