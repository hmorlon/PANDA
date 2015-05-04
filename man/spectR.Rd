\name{spectR}
\alias{spectR}

\title{
Spectral density plot of a phylogeny
}

\description{
Computes the spectra of eigenvalues for the modified graph Laplacian of a phylogenetic tree, identifies the spectral gap, then convolves the eigenvalues with a Gaussian kernel, and plots them alongside all eigenvalues ranked in descending order.}

\usage{
spectR(phylo, method=c("standard","normal"))
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
  \item{method}{
whether the spectral density plot will be normalized or not
}

\options{
the "normal" method normalizes the graph Laplacian by the diagonal matrix.	

\value{
	Density plot of the spectra of eigenvalues computed for a phylogeny.
}

\author{E Lewitus}


\examples{
require(geiger)
require(igraph)	
data(Cetacea)
spectR(Cetacea,method="standard")
}
