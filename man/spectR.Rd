\name{spectR}
\alias{spectR}

\title{
Spectral density plot of a phylogeny
}

\description{
Computes the spectra of eigenvalues for the modified graph Laplacian of a phylogenetic tree, identifies the spectral gap, then convolves the eigenvalues with a Gaussian kernel, and plots them alongside all eigenvalues ranked in descending order.
}

\usage{
spectR(phylo, method=c("standard"))
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{method}{
the method used to compute the spectral density, which can either be "standard", "normal1", or "normal2". if set to "normal1", compute the spectral density normalized to the degree matrix. if set to "normal2", compute the spectral density normalized to the number of eigenvalues. if set to "standard", compute the unnormalized version of the spectral density (see the associated paper for an explanation)}
}

\value{
	a list with the following components
	\item{eigenvalues} {the vector of eigenvalues} 
	\item{eigengap} {the position of the largest difference between eigenvalues, giving the number of modalities in the tree} 
	}

\details{
The function also returns a plot of the spectral density and a plot of all eigenvalues ranked in descending order.
}

\author{E Lewitus}

\references{
	Lewitus, E., Morlon, H., A graph-theoretical approach to characterising and comparing phylogenies (submitted) }


\examples{
data(Cetacea)
spectR(Cetacea,method="standard")
}
