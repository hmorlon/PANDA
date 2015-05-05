\name{JSDtree}
\alias{JSDtree}

\title{
The Jensen-Shannon distance metric between phylogenies
}

\description{
Computes the Jensen-Shannon distance metric between spectral density profiles of phylogenies.
}

\usage{
JSDtree(trees,method=c("standard"))
}

\arguments{
  \item{phylo}{
a list of objects of type 'phylo' (see ape documentation)}
   \item{method}{
whether the spectral density will be normalized or not}
}

\details{
	It is possible to specify whether a standard profile is returned, a profile normalized to the degree matrix ("normal1"), or a profile normalized to the number of eigenvalues ("normal2").
}

\value{
	Writes the Jensen-Shannon distance metrics between phylogeny pairs to a text file and plots a hierarchical cluster and heatmap based on the distance metrics calculated between trees.
}

\author{E Lewitus}


\examples{
require(TESS)
require(igraph)	
sim.globalBiDe.age(n=10,age=10,0.15,0.05,MRCA=TRUE) -> trees
JSDtree(trees)
}
