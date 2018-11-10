
\name{true.div}
\alias{true.div}
\title{True diversity}
\description{Compute true diversities of individual samples from a vector of relative abundances or a matrix (OTU table) containing multiple samples.}
\usage{true.div(abund,qvalue,tree)}
\arguments{
  \item{abund}{A vector or a matrix indicating the relative abundances of one or multiple samples. If a matrix is provided, columns should be samples and rows OTUs. The relative abundances of each sample (vector or each column of the matrix) must sum to 1.}
  \item{qvalue}{A positive number (>0). It can be an integer or contain decimals.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
}
