\name{trait_spectR}
\alias{trait_spectR}

\title{
Spectral density plot of phylogeny trait data
}

\description{
Computes the spectra of eigenvalues for the modified graph Laplacian of a phylogenetic tree with associated tip data, convolves the eigenvalues with a Gaussian kernel and plots the density profile of eigenvalues, and estimates the summary statistics of the profile.}

\usage{
trait_spectR(phylo, data, plot=F)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{data}{
a vector of trait data associated with the tips of the phylo object; tips and trait data should be aligned}
}

\value{
	a list with the following components:
	\item{eigenvalues}{
the vector of eigenvalues} 
	\item{splitter}{
the largest (or principal) eigenvalue of the spectral density profile}
\item{fragmenter}{
the skewness of the spectral density profile} 
	\item{tracer}{
the largest y-axis value of the spectral density profile}
	}



\references{
	Lewitus, E., Morlon, H., Characterizing and comparing phylogenies from their Laplacian spectrum, bioRxiv doi: http://dx.doi.org/10.1101/026476}

\author{E Lewitus}


\examples{
tr<-rtree(10)
dat<-runif(10,1,2)
trait_spectR(tr,dat,plot=T)
}
