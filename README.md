# hilldiv 1.4 (Last update: August 2019)

**hilldiv** is an R package that provides a set of functions to assist analysis of diversity for diet reconstruction, microbial community profiling or more general ecosystem characterisation analyses based on Hill numbers, using OTU tables and associated phylogenetic trees as inputs. The package includes functions for (phylo)diversity measurement, (phylo)diversity profile plotting, (phylo)diversity comparison between samples and groups,  (phylo)diversity partitioning and (dis)similarity measurement. All of these grounded in abundance-based and incidence-based Hill numbers.

The statistical framework developed around Hill numbers encompasses many of the most broadly employed diversity (e.g. richness, Shannon index, Simpson index), phylogenetic diversity (e.g. Faith’s PD, Allen’s H, Rao’s quadratic entropy) and dissimilarity (e.g. Sørensen index, Unifrac distances) metrics. This enables the most common analyses of diversity to be performed while grounded in a single statistical framework. For details about the use of Hill numbers in molecularly characterised biological systems, read the following article.

Alberdi A, Gilbert MTP. (2019). A guide to the application of Hill numbers to DNA‐based diversity analyses. *Molecular Ecology Resources*. 19(4): 804-817. [https://doi.org/10.1111/1755-0998.13014](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014)

**If you want/need to cite hilldiv:** Alberdi A, Gilbert MTP. 2019. hilldiv: an R package for the integral analysis of diversity based on Hill numbers. bioRxiv 545665. https://www.biorxiv.org/content/10.1101/545665v1

Some recent papers that used hilldiv:

* Alberdi et al. (2019) Diet diversification shapes broad-scale distribution patterns in European bats. *bioRxiv*. https://doi.org/10.1101/704759

* Siren et al. (2019) Taxonomic and Functional Characterization of the Microbial Community During Spontaneous in vitro Fermentation of Riesling Must. *Frontiers in Microbiology*. https://doi.org/10.3389/fmicb.2019.00697

* Chang et al. (2019) The roles of morphological traits, resource variation and resource partitioning associated with the dietary niche expansion in the fish‐eating bat *Myotis pilosus*. *Molecular Ecology*. https://doi.org/10.1111/mec.15127

# List of functions
| Function | Short explanation | Documentation |
| ------------- | ------------- | ------------- |
| hill.div()  | Neutral or phylogenetic Hill number computation   | [LINK](documentation/hill.div.md) |
| index.div()  | Neutral or phylogenetic diversity index computation   | [LINK](documentation/index.div.md) |
| div.profile()  | (Phylo)Diversity profiles of individual samples or groups of samples | [LINK](documentation/div.profile.md) |
| div.test()  | Diversity comparison between two or multiple groups of samples   | [LINK](documentation/div.test.md) |
| div.test.plot()  | Visual representation of div.test() | [LINK](documentation/div.test.plot.md) |
| depth.cov()  | Assessment of the sequencing depth per sample   | [LINK](documentation/depth.cov.md) |
| div.part()  | Hierarchical diversity partitioning   | [LINK](documentation/div.part.md) |
| alpha.div()  | Alpha diversity computation   | [LINK](documentation/alpha.div.md) |
| gamma.div()  | Gamma diversity computation   | [LINK](documentation/gamma.div.md) |
| beta.dis()  | (Dis)similarity computation based on beta diversities   | [LINK](documentation/beta.dis.md) |
| pair.dis()  | Pairwise (dis)similarity computation based on beta diversities   | [LINK](documentation/pair.dis.md) |
| pair.dis.plot()  | Visual representation of pair.dis()   | [LINK](documentation/pair.dis.plot.md) |
| UqN()  | Jaccard-type overlap computation from beta diversities | [LINK](documentation/UqN.md) |
| CqN()  | Sørensen-type overlap from beta diversities   | [LINK](documentation/CqN.md) |
| SqN()  | Jaccard-type turnover-complement from beta diversities   | [LINK](documentation/SqN.md) |
| VqN()  | Sørensen-type turnover-complement from beta diversities   | [LINK](documentation/VqN.md) |
| match.data()  | Filter OTU tables and trees to match OTUs in both data files   | [LINK](documentation/match.data.md) |
| depth.filt()  | Filter samples according to a minimum sequencing depth threshold   | [LINK](documentation/depth.filt.md) |
| copy.filt()  | Filter OTUs according to a minimum copy number threshold   | [LINK](documentation/copy.filt.md) |
| tss()  | Total sum scaling per sample   | [LINK](documentation/tss.md) |
| tree.depth()  | Tree depth computation   | [LINK](documentation/tree.depth.md) |

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

### v1.4 | August 2019
- Implementation of roxygen2 

### v1.2.3 | July 2019
- Added the function index.div() to compute diversity indices related to Hill numbers.
- Added the match.data() function.
- Added the depth.filt() function.
- Added the auxiliary tree.depth() function.
- Added sample weighting option and corrected "object 'weight.L2' not found" error in pair.dis() function.
- Added automatic tss normalisation to alpha.div() and gamma.div() functions.

### v1.2.2 | June 2019
- Added copy.filt() function to filter OTUs according to absolute or relative copy number thresholds.

### v1.2.1 | May 2019
- Option to plot pairwise mean comparison statistical significance values added to pair.dist.plot().
- Added depth.cov() function for assessment of the sequencing depth per sample based on observed and estimated Hill numbers.

### v1.1 | March 2019
- Dependent "geiger" package added to the Description document.
- Magnify issue corrected in pair.dis.plot.r

# Workflow
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/hilldiv_workflow.png" title="Hilldiv workflow">

# Documentation
I am creating extensive documentation about Hill numbers and implementation of *hilldiv* in the [Hilldiv WIKI](https://github.com/anttonalberdi/hilldiv/wiki)

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
names(otu.vector) <- rownames(otu.table)
hierarchy.table <- bat.diet.hierarchy
tree <- bat.diet.tree
````
## Preliminary analyses and filtering
- copy.filt()
- depth.cov()
- match.data()

## Diversity measurement and visualisation
- hill.div()
- index.div()
- div.profile()

## Diversity comparison
- div.test()
- div.test.plot()

## Diversity partitioning
- div.part()
- alpha.div()
- gamma.div()

## Dissimilarity measurement and visualisation
- beta.dis()
- pair.dis()
- pair.dis.plot()
- CqN()
- UqN()
- SqN()
- VqN()

## Auxiliary functions
- tss()
- tree.depth()

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
