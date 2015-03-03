\name{JSDtree}
\alias{JSDtree}

\title{
The Jensen-Shannon divergence metric between phylogenies
}

\description{
Compute the square-root of the Jensen-Shannon divergence metric between spectral density profiles of phylogenies.
}

\usage{
spectR(phylo)
}

\arguments{
  \item{phylo}{
a list of objects of type 'phylo' (see ape documentation)
   \item{method}{
whether the spectral density will be normalized or not
}


\value{
	Writes the square-root of Jensen-Shannon divergence metrics between phylogeny pairs to a text file and plots a hierarchical cluster and heatmap based on the divergence metrics calculated between trees.
}

\author{E Lewitus}


\examples{
require(TESS)
require(igraph)	
sim.globalBiDe.age(n=10,age=10,0.15,0.05,MRCA=TRUE) -> trees
JSDtree(trees)
}
