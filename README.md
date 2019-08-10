# hilldiv 1.2.3 (Last update: July 2019)

**hilldiv** is an R package that provides a set of functions to assist analysis of diversity for diet reconstruction, microbial community profiling or more general ecosystem characterisation analyses based on Hill numbers, using OTU tables and associated phylogenetic trees as inputs. The package includes functions for (phylo)diversity measurement, (phylo)diversity profile plotting, (phylo)diversity comparison between samples and groups,  (phylo)diversity partitioning and (dis)similarity measurement. All of these grounded in abundance-based and incidence-based Hill numbers.

The statistical framework developed around Hill numbers encompasses many of the most broadly employed diversity (e.g. richness, Shannon index, Simpson index), phylogenetic diversity (e.g. Faith’s PD, Allen’s H, Rao’s quadratic entropy) and dissimilarity (e.g. Sørensen index, Unifrac distances) metrics. This enables the most common analyses of diversity to be performed while grounded in a single statistical framework. For details about the use of Hill numbers in molecularly characterised biological systems, read the following article.

Alberdi A, Gilbert MTP. (2019). A guide to the application of Hill numbers to DNA‐based diversity analyses. *Molecular Ecology Resources*. 19(4): 804-817. [https://doi.org/10.1111/1755-0998.13014](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014)

**If you want/need to cite hilldiv:** Alberdi A, Gilbert MTP. 2019. hilldiv: an R package for the integral analysis of diversity based on Hill numbers. bioRxiv 545665. https://www.biorxiv.org/content/10.1101/545665v1

Some recent papers that used hilldiv:

* Alberdi et al. (2019) Diet diversification shapes broad-scale distribution patterns in European bats. *bioRxiv*. https://doi.org/10.1101/704759 

* Siren et al. (2019) Taxonomic and Functional Characterization of the Microbial Community During Spontaneous in vitro Fermentation of Riesling Must. *Frontiers in Microbiology*. https://doi.org/10.3389/fmicb.2019.00697

* Chang et al. (2019) The roles of morphological traits, resource variation and resource partitioning associated with the dietary niche expansion in the fish‐eating bat *Myotis pilosus*. *Molecular Ecology*. https://doi.org/10.1111/mec.15127

# List of functions
| Function | Short explanation |
| ------------- | ------------- |
| hill.div()  | Neutral or phylogenetic Hill number computation   |
| index.div()  | Neutral or phylogenetic diversity index computation   |
| div.profile()  | (Phylo)Diversity profiles of individual samples or groups of samples   |
| div.test()  | Diversity comparison between two or multiple groups of samples   |
| div.test.plot()  | Visual representation of div.test() |
| depth.cov()  | Assessment of the sequencing depth per sample   |
| div.part()  | Hierarchical diversity partitioning   |
| alpha.div()  | Alpha diversity computation   |
| gamma.div()  | Gamma diversity computation   |
| beta.dis()  | (Dis)similarity computation based on beta diversities   |
| pair.dis()  | Pairwise (dis)similarity computation based on beta diversities   |
| pair.dis.plot()  | Visual representation of pair.dis()   |
| UqN()  | Jaccard-type overlap computation from beta diversities |
| CqN()  | Sørensen-type overlap from beta diversities   |
| SqN()  | Jaccard-type turnover-complement from beta diversities   |
| VqN()  | Sørensen-type turnover-complement from beta diversities   |
| copy.filt()  | Filter OTUs according to a minimum copy number threshold   |
| tss()  | Total sum scaling per sample   |
| tree.depth()  | Tree depth computation   |

# Installation
To install **hilldiv** in your R environment, you need to 1) install devtools, 2) load devtools library, 3) install **hilldiv** using devtools and 4) finally load **hilldiv** library to your environment.

````R
install.packages("devtools")
library(devtools)
install_github("anttonalberdi/hilldiv")
library(hilldiv)
````
If not installed, it will automatically install the following dependencies: ggplot2, ggpubr, RColorBrewer, data.table, ape, ade4, iNEXT, iNextPD.

If the console returns the following error:
````R
"tar: Failed to set default locale"
````
Type the following in the console and restart R
````R
system('defaults write org.R-project.R force.LANG en_US.UTF-8')
````

# Changelog

### v1.2.3 | July 2019
- Added the function index.div() to compute diversity indices related to Hill numbers.
- Added the depth.filt() function.
- Added the auxiliary tree.depth() function.
- Added sample weighting option to pair.dis() function.
- Added automatic tss normalisation to alpha.div() and gamma.div() functions.

### v1.2.2 | June 2019
- Added copy.filt() function to filter OTUs according to absolute or relative copy number thresholds.

### v1.2.1 | May 2019
- Option to plot pairwise mean comparison statistical significance values added to pair.dist.plot().
- Added depth.cov() function for assessment of the sequencing depth per sample based on observed and estimated Hill numbers.

### v1.1 | March 2019
- Dependent "geiger" package added to the Description document.
- Magnify issue corrected in pair.dis.plot.r

# Applications and functions
Note that detailed information about the use of *hilldiv* can be found in the [Hilldiv WIKI](https://github.com/anttonalberdi/hilldiv/wiki)
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
names(otu.vector) <- rownames(otu.table)
hierarchy.table <- bat.diet.hierarchy
tree <- bat.diet.tree
````
## Preliminary analyses and filtering
Before begining the diversity analyses it is recommendable to assess the reliability and representativeness of the data and filter OTUs to get rid of potentially artifactual sequences. **hilldiv** includes a few functions to with that purpose.
### copy.filt()
It is highly probable that OTUs with very low representation (number of DNA sequence copies) are artifactual, rather than real biological sequences. Therefore, it is recommendable to remove OTUs represented by a number of sequences below a certain threshold. Although the use of absolute thresholds (e.g. remove all sequences with less than 10 copies) is the most commonly employed strategy, ideally thresholds relative to the sequencing depth should be employed (Alberdi et al. 2018) in order to apply a comparable treatment to samples characterised with a different sequencing depth. The function copy.filt() enables applying both absolute (if threshold value is an integer) and relative (if threshold value is a decimal number) copy number thresholds

````R
#Remove singletons from one sample (absolute threshold)
otu.vector.filtered <- copy.filt(otu.vector,2)
#Remove singletons from an OTU table (absolute threshold)
otu.table.filtered <- copy.filt(otu.table,2)
#Remove OTUs with less than 0.1% of the total number of sequences (sequencing depth) per sample (relative threshold)
otu.table.filtered <- copy.filt(otu.table,0.001)
````
### depth.cov()
Assessing whether the sequencing depth of each sample is enough to recover the entire diversity of a sample is an important step to ensure reliable and unbiased comparisons across samples. The function depth.cov() relies on diversity estimations based on Hill numbers (Chao & Jost 2015) to calculate the percentage of estimated diversity covered in each sample.

````R
#Depth coverage assessment of multiple samples based on the order of diversity 0
depth.cov(otu.table,0)
#Depth coverage assessment of a single sample based on the order of diversity 1
depth.cov(otu.vector,1)
````

## Diversity measurement and visualisation
### index.div()
Neutral (richness, Shannon index, Simpson index) or phylogenetic (Faith's PD, Allen's H, Rao's Q) diversity indices related to Hill numbers from a vector object (one sample) or an OTU table (matrix or data.frame object; multiple samples). A ultrametric tree object (phylo) is necessary to compute phylogenetic diversity indices. Note that if using a tree the tip labels and the 'names' (vectors) or 'rownames' (matrices) need to be identical. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. If the vector or the OTU table columns do not sum to 1, the data is TSS-normalised.
````R
index.div(otu.vector)
index.div(otu.table,index="shannon")
index.div(otu.table,tree=tree)
index.div(otu.table,tree=tree,index="allen")
index.div(otu.table,tree,"rao")
````

### hill.div()
Neutral or phylogenetic Hill numbers computation from a vector object (one sample) or an OTU table (matrix or data.frame object; multiple samples). Providing the tree argument yields phylodiversity values. Note that if using a tree the tip labels and the 'names' (vectors) or 'rownames' (matrices) need to be identical. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. If the vector or the OTU table columns do not sum to 1, the data is TSS-normalised.
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
Visual comparison between the diversity levels of two or multiple groups of samples. The 'chart' argument enables selecting between boxplot, violin plot and jitter plot. The 'stat' argument enables plotting pairwise mean comparisons across groups.
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

#Pairwise mean comparison statistical significances
div.test.plot(contrast.div.q0,chart="jitter",stat=TRUE,flip=TRUE)
div.test.plot(contrast.div.q0,stat=TRUE,comb=list(c("Myotis myotis","Myotis capaccinii")),symbol=TRUE)
````
<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.test.plot.png" width="400" title="div.test.plot() with pairwise comparisons">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.test.plot.violin.png" width="400" title="Violin div.test.plot() with pairwise comparisons">

## Sequencing depth assessment
### depth.cov()
Assessment of the sequencing depth per sample. Computes both observed and estimated (based on Chao & Jost 2015) Hill numbers, and calculates the percentage of estimated diversity coverage per sample. The function requires absolute abundances (rather than relative) to be inputed, and estimations are most accurate in the range of 0<q<3.

````R
depth.cov(otu.table,0)
depth.cov(otu.vector,2)
depth.cov(otu.table,qvalue=1)
````

## Diversity partitioning
### div.part
Hierarchical diversity partitioning. The function assumes a 2-level hierarchy and yields alpha, gamma and beta values based on abundance data. If a hierarchy table is provided, the function yields alpha, gamma and beta values based on incidence data, i.e. alpha diversity reflects the incidence-based diversity of groups.
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
For functions hill.div, alpha.div, gamma.div, div.part...

|       | Sample1 | Sample2 |Sample3  |SampleN  |
| ------------- | ------------- | ------------- | ------------- |------------- |
| OTU1  | 604   |543    |204    |...    |
| OTU2  | 50   |210    |450    |...  |
| OTU3  | 153   |105    |314    |...  |
| OTUN  | ...   |...    |...    |...  |

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
* Alberdi A., Gilbert M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19(4), 804-817.
* Chao, A. & Jost, L. (2015) Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution, 6, 873-882.
* Chao, A., Chiu, C.-H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93(9), 2037–2051.
* Chao, A., Chiu, C.-H., & Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences, 365(1558), 3599–3609.
* Hill, M. O. (1973). Diversity and Evenness: A Unifying Notation and Its Consequences. Ecology, 54(2), 427–432.
* Hsieh, T. C., & Chao, A. (2017). Rarefaction and Extrapolation: Making Fair Comparison of Abundance-Sensitive Phylogenetic Diversity among Multiple Assemblages. Systematic Biology, 66(1), 100–111.
* Hsieh, T. C., Ma, K. H., & Chao, A. (2016). iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution / British Ecological Society, 7(12), 1451–1456.
* Jost, L. (2006). Entropy and diversity. Oikos , 113, 363–375.
* Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88(10), 2427–2439.
* Marcon, E., & Hérault, B. (2015). entropart: An R Package to Measure and Partition Diversity. Journal of Statistical Software, Articles, 67(8), 1–26.
* Tuomisto, H. (2010). A diversity of beta diversities: straightening up a concept gone awry. Part 1. Defining beta diversity as a function of alpha and gamma diversity. Ecography, 33(1), 2–22.
