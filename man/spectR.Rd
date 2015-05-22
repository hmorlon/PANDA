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
the method used to compute the spectral density. "standard" corresponds to the spectral density with no normalization. "normal1" corresponds to a normalization by the degree matrix. "normal2" corresponds to a normalization by the number of eigenvalues.}
}

\value{
	a list with the following components
	\item{eigenvalues} {the vector of eigenvalues} 
	\item{eigengap} {the position of the largest difference between eigenvalues} 
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
