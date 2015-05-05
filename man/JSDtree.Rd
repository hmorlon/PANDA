\name{JSDtree}
\alias{JSDtree}

\title{
The Jensen-Shannon distance metric between phylogenies
}

\description{
Computes the Jensen-Shannon distance metric between spectral density profiles of phylogenies.
}

\usage{
JSDtree(phylo,method=c("standard","normal"))
}

\arguments{
  \item{phylo}{
a list of objects of type 'phylo' (see ape documentation)}
   \item{method}{
whether the spectral density will be normalized or not}
}


\value{
	Writes the Jensen-Shannon distance metrics between phylogeny pairs to a text file and plots a hierarchical cluster and heatmap based on the distance metrics calculated between trees.
}

\author{E Lewitus}


\examples{
sim.globalBiDe.age(n=10,age=10,0.15,0.05,MRCA=TRUE) -> trees
JSDtree(trees)
}
