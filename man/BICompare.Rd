\name{BICompare}
\alias{BICompare}

\title{
Compute the confidence of a modality in a phylogeny
}

\description{
Computes the BIC values for a specified number of modalities in the distance matrix of a tree and a control matrix. 
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

\details{
	Note that the "poisson" option should only be used in certain cases (e.g., when the tree is extremely imbalanced)}


\author{E Lewitus}


\examples{
data(Cetacea)
BICompare(Cetacea,5,method="gaussian")
}
