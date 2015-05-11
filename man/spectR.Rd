\name{spectR}
\alias{spectR}

\title{
Spectral density plot of a phylogeny
}

\description{
Computes the spectra of eigenvalues for the modified graph Laplacian of a phylogenetic tree, identifies the spectral gap, then convolves the eigenvalues with a Gaussian kernel, and plots them alongside all eigenvalues ranked in descending order.
}

\usage{
spectR(phylo, method=c("standard","normal1","normal2"))
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{method}{
whether the spectral density plot will be normalized or not}
}

\details{
	It is possible to specify whether a standard profile is returned, a profile normalized to the degree matrix ("normal1"), or a profile normalized to the number of eigenvalues ("normal2").
}

\value{
	List of eigenvalues
	Eigengap
	Density plot of the spectra of eigenvalues computed for a phylogeny.
}

\author{E Lewitus}


\examples{
data(Cetacea)
spectR(Cetacea,method="standard")
}
