\name{BICompare}
\alias{BICompare}

\title{
Compute the confidence of modalities in a phylogeny
}

\description{
Computes the BIC values for different numbers of modalities in the distance matrix of a tree and a control matrix. 
}

\usage{
BICompare(phylo, t, method=c("gaussian","poisson"))
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{t}{
the number of modalities to be tested}
  \item{method}{
whether the control matrix will be drawn from a poisson or gaussian distribution}	
}


\value{
	BIC values for finding t modalities in the distance matrix of a tree and a control matrix.
}

\author{E Lewitus}


\examples{
data(Cetacea)
BICompare(Cetacea,5,method="gaussian")
}
