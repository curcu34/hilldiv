**hilldiv** is an R package that provides a set of functions to assist analysis of diversity for diet reconstruction, microbial community profiling or more general ecosystem characterisation analyses based on Hill numbers, using OTU tables and associated phylogenetic trees as inputs. The package includes functions for (phylo)diversity measurement, (phylo)diversity profile plotting, (phylo)diversity comparison between samples and groups,  (phylo)diversity partitioning and (dis)similarity measurement. All of these grounded in abundance-based and incidence-based Hill numbers.

The statistical framework developed around Hill numbers encompasses many of the most broadly employed diversity (e.g. richness, Shannon index, Simpson index), phylogenetic diversity (e.g. Faith’s PD, Allen’s H, Rao’s quadratic entropy) and dissimilarity (e.g. Sørensen index, Unifrac distances) metrics. This enables the most common analyses of diversity to be performed while grounded in a single statistical framework. 

* Installation
* Applications and functions
* References

# Installation
To install **hilldiv** in your R environment, you need to 1) install devtools, 2) load devtools library, 3) install **hilldiv** using devtools and 4) finally load **hilldiv** library to your environment.

````R
install.packages("devtools")
library(devtools)
install_github("anttonalberdi/hilldiv")
library(hilldiv)
````
If not installed, it will automatically install the following dependencies: ggplot2, RColorBrewer, data.table, ape, ade4, iNEXT, iNextPD.

# Applications and functions
## Data
All the different applications and functions are reproduced with the data included in this package. 
````R
#Load data
data(bat.diet.otutable)
data(bat.diet.hierarchy)
data(bat.diet.tree)

#Create simple objects
otu.table <- bat.diet.otutable
otu.vector <- bat.diet.otutable[,1]
hierarchy.table <- bat.diet.hierarchy
tree <- bat.diet.tree
````
## Diversity measurement and visualisation
### hill.div()
Neutral or phylogenetic Hill number computation from a vector object (one sample) or an OTU table (matrix or data.frame object; multiple samples). Providing the tree argument yields phylodiversity values. Note that if using a tree the tip labels and the 'names' (vectors) or 'rownames' (matrices) need to be identical. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. If the vector or the OTU table columns do not sum to 1, the data is TSS-normalised.
````R
hill.div(otu.vector,0)
hill.div(otu.table,1)
hill.div(otu.table,qvalue=2)
hill.div(otu.table,1,tree)
hill.div(otu.table,1,tree,type="incidence")
hill.div(otu.table,qvalue=2,tree=tree)
````
### div.profile() - chart
(Phylo)Diversity profiles of individual samples or groups of samples. Diversity profiles show the relation between the order of diversity (q-value) and the respective Hill numbers, thus providing information about the richness and evenness of a sample at a glance.

````R
# One sample
div.profile(otu.vector)

# Multiple individual samples (first 5 samples of the OTU table)
div.profile(otu.table[,c(1:5)])

# Multiple groups (aggregated samples)
div.profile(otu.table,hierarchy=hierarchy.table,colour=c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7"))
````
<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.one.png" width="350" title="One sample">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.multiple.png" width="350" title="Multiple samples">

## Diversity comparison
### div.test()
Diversity comparison between two or multiple groups of samples. If the tree argument is used the test compares phylodiversity values. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. 

````R
#USAGE#
#Contrast based on Hill numbers
div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
#Contrast based on phylogenetic Hill numbers
div.test(otu.table,qvalue=1,hierarchy=hierarchy.table,tree=tree)
````
### div.test.plot() - chart
Visual comparison between the diversity levels of two or multiple groups of samples. The chart argument enables selecting between boxplot and jitter plot.
````R
contrast.div.q0 <- div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
colours <- c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7")

#Box plot
div.test.plot(contrast.div.q0)
div.test.plot(contrast.div.q0,chart="box")

#Jitter plot
div.test.plot(contrast.div.q0,chart="jitter",colour=colours)

#Violin plot
div.test.plot(contrast.div.q0,chart="violin",colour=c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7"))
````

## Diversity partitioning
### div.part
The function assumes a 2-level hierarchy and yields alpha, gamma and beta values based on abundance data. If a hierarchy table is provided, the function yields alpha, gamma and beta values based on incidence data, i.e. alpha diversity reflects the incidence-based diversity of groups.
````R
#Abundance-based
div.part(otu.table,qvalue=1)
div.part(otu.table,qvalue=1,type="abundance")
div.part(otu.table,qvalue=0,tree=tree)

#Incidence-based
div.part(otu.table,qvalue=0,type="incidence",hierarchy=hierarchy.table)
````

## Dissimilarity measurement and visualisation
###  beta.dis()
The function beta.dis() performs similarity or dissimilarity measurement based on Hill numbers beta diversity, sample size and order of diversity. The function can be run by inputting those values manually, or by using the list object outputted by the div.part() function, which contains all the mentioned information. As specified by the argument “metric”: the function can compute the following similarity measures: the Sørensen-type overlap (CqN), the Jaccard-type overlap (UqN), the Sørensen-type turnover-complement (VqN), and the Jaccard-type turnover-complement (SqN). The argument ‘type’ enables either similarities or dissimilarities (one-complements of the similarity values) to be outputted.
````R
#Using custom, beta-diversity, q-value and samples sizes
beta.dis(beta=4.5,qvalue=1,N=8)
beta.dis(beta=4.5,qvalue=1,N=8,metric="C",type="similarity")

#Computing from an div.part() derived object
div.part.object  <- div.part(otu.table,qvalue=0,tree=tree)
beta.dis(div.part.object)
beta.dis(div.part.object,metric="S",type="similarity")
````
###  pair.dis()
The function pair.dis() performs pairwise diversity partitioning and yields matrices containing pairwise beta diversity and (dis)similarity measures. If a hierarchy table is provided, pairwise calculations can be carried out at some or all the specified hierarchical levels. The results are outputted as a list of matrices.
````R
pair.dis(otu.table,qvalue=0,hierarchy=hierarchy.table)
pair.dis(otu.table,qvalue=0,hierarchy=hierarchy.table,level="2")
````

###  pair.dis.plot() - chart
The related function pair.dis.plot() uses any of the dissimilarity matrices yielded by pair.dis() (e.g. 1-UqN) to visualize it either as a NMDS chart, a qgraph plot or a heatmap/correlogram.

````R
pair.div.q0.L2 <- pair.dis(otu.table[,sort(colnames(otu.table))],qvalue=0,hierarchy=hierarchy.table,level="2")
#Plot NMDS
pair.dis.plot(pair.div.q0.L2$L2_CqN,hierarchy=hierarchy.table,type="NMDS",level=2)
#Plot qgraph
pair.dis.plot(pair.div.q0.L2$L2_CqN,hierarchy=hierarchy.table,type="qgraph",level=2,magnify=TRUE)
````

## Auxiliary functions
###  tss()
Performs total sum scaling.
````R
tss(otu.table)
tss(otu.vector)
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

| Sample | Group |
| ------------- | ------------- |
| Sample1  | Group1   |
| Sample2  | Group1   |
| Sample3  | Group2   |
| Sample4  | Group2   |
| Sample5  | Group3   |
| ...  | ...   |
| SampleN  | GroupN   |

* Sample names and OTU table row names must match
* Groups and supergroups must be nested (all samples in a group must be included in the same supergroup)

# References
* Chao, A., Chiu, C.-H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93(9), 2037–2051.
* Chao, A., Chiu, C.-H., & Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences, 365(1558), 3599–3609.
* Hill, M. O. (1973). Diversity and Evenness: A Unifying Notation and Its Consequences. Ecology, 54(2), 427–432.
* Hsieh, T. C., & Chao, A. (2017). Rarefaction and Extrapolation: Making Fair Comparison of Abundance-Sensitive Phylogenetic Diversity among Multiple Assemblages. Systematic Biology, 66(1), 100–111.
* Hsieh, T. C., Ma, K. H., & Chao, A. (2016). iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution / British Ecological Society, 7(12), 1451–1456.
* Jost, L. (2006). Entropy and diversity. Oikos , 113, 363–375.
* Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88(10), 2427–2439.
* Marcon, E., & Hérault, B. (2015). entropart: An R Package to Measure and Partition Diversity. Journal of Statistical Software, Articles, 67(8), 1–26.
* Tuomisto, H. (2010). A diversity of beta diversities: straightening up a concept gone awry. Part 1. Defining beta diversity as a function of alpha and gamma diversity. Ecography, 33(1), 2–22.
