\name{gamma.curve}
\alias{gamma.curve}
\title{Gamma curve}
\usage{}
\arguments{
  \item{otutable}{A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs. The relative abundances of each sample (columns) must sum to 1.}
  \item{qvalue}{A positive number (>0). It can be an integer or contain decimals.}
  \item{tree}{A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).}
  \item{hierarchy}{
  \item{steps}{ 
  \item{iter}{ Integer | The number of iterations to create the gamma curve. The higher the number, the lower is the standard error and the smoother is the curve. iter = 50}
  \item{summary}{ Logical (TRUE/FALSE) | If 'TRUE' the summary data used to generate the gamma curve is also printed. }

  }
