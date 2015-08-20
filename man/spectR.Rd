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
	a list with the following components:
	\item{eigenvalues}{
the vector of eigenvalues} 
	\item{principal_eigenvalue}{
the largest (or principal) eigenvalue of the spectral density profile}
\item{asymmetry}{
the skewness of the spectral density profile} 
	\item{peakedness}{
the kurtosis of the spectral density profile}
	\item{eigengap}{
the position of the largest difference between eigenvalues, giving the number of modalities in the tree}
	}
	
\details{
	Note that the eigengap should in principle be computed with the "standard" option}


\references{
	Lewitus, E., Morlon, H., A graph-theoretical approach to characterising and comparing phylogenies (submitted) }

\author{E Lewitus}

\seealso{
\code{\link{plot_spectR}}, \code{\link{JSDtree}}, \code{\link{BICompare}}
}


\examples{
data(Cetacea)
spectR(Cetacea,method="standard")
}