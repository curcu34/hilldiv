# index.div
Neutral (richness, Shannon index, Simpson index) or phylogenetic (Faith's PD, Allen's H, Rao's Q) diversity indices related to Hill numbers from a vector object (one sample) or an OTU table (matrix or data.frame object; multiple samples). A ultrametric tree object (phylo) is necessary to compute phylogenetic diversity indices. Note that if using a tree the tip labels and the 'names' (vectors) or 'rownames' (matrices) need to be identical. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. If the vector or the OTU table columns do not sum to 1, the data is TSS-normalised.
````R
index.div(otu.vector)
index.div(otu.table,index="shannon")
index.div(otu.table,tree=tree)
index.div(otu.table,tree=tree,index="allen")
index.div(otu.table,tree,"rao")
````
