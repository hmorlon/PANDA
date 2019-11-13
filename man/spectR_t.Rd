\name{spectR_t}
\alias{spectR_t}

\title{
Spectral density plot of phylogenetic trait data
}

\description{
Computes the spectra of eigenvalues for the modified graph Laplacian of a phylogenetic tree with associated tip data, convolves the eigenvalues with a Gaussian kernel and plots the density profile of eigenvalues, and estimates the summary statistics of the profile.}

\usage{
spectR_t(phylo, dat, draw=F)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{dat}{
a vector of trait data associated with the tips of the phylo object; tips and trait data should be aligned}
   \item{draw}{
   	if true, the spectral density profile of the phylogenetic trait data is plotted
   }	
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
	Lewitus, E., Morlon, H. (2019) Characterizing and comparing phylogenetic trait data from their normalized Laplacian spectrum, bioRxiv doi: https://doi.org/10.1101/654087}


\author{E Lewitus}


\examples{
tr<-rtree(10)
dat<-runif(10,1,2)
spectR_t(tr,dat,draw=TRUE)
}
