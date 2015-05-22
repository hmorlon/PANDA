\name{JSDtree}
\alias{JSDtree}

\title{
The Jensen-Shannon distance metric between phylogenies
}

\description{
Computes the Jensen-Shannon distance metric between spectral density profiles of phylogenies.
}

\usage{
JSDtree(phylo,meth=c("standard"))
}

\arguments{
  \item{phylo}{
a list of objects of type 'phylo' (see ape documentation)}
   \item{meth}{
the method used to compute the spectral density, which can either be "standard", "normal1", or "normal2". if set to "normal1", compute the spectral density normalized to the degree matrix. if set to "normal2", compute the spectral density normalized to the number of eigenvalues. if set to "standard", compute the unnormalized version of the spectral density (see the associated paper for an explanation)}
}


\value{
	A matrix providing the Jensen-Shannon distance values between phylogeny pairs. 
	The function also plots a hierarchical cluster and heatmap based on the distance metrics calculated between trees.
}


\references{
Lewitus E., Morlon H. (submitted) A graph-theoretical approach to characterising and comparing phylogenies
}

\author{E Lewitus}


\examples{
trees<-sim.globalBiDe.age(n=10,age=10,0.15,0.05,MRCA=TRUE)
JSDtree(trees)
}
