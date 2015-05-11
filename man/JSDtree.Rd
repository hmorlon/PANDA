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
the method to be used. either "standard" or "normal". if set to "normal", compute the normalized version of the spectral density (see the paper for an explanation)}
}

\details{
	It is possible to specify whether a standard profile is returned, a profile normalized to the degree matrix ("normal1"), or a profile normalized to the number of eigenvalues ("normal2").
}

\value{
	Writes the Jensen-Shannon distance metrics between phylogeny pairs to a text file and plots a hierarchical cluster and heatmap based on the distance metrics calculated between trees.
}

\author{E Lewitus}


\examples{
sim.globalBiDe.age(n=10,age=10,0.15,0.05,MRCA=TRUE) -> trees
JSDtree(trees)
}
