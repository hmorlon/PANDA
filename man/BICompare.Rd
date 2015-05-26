\name{BICompare}
\alias{BICompare}

\title{
Compute the confidence of a modality in a phylogeny
}

\description{
Computes the BIC values for a specified number of modalities in the distance matrix of a tree and a control matrix. 
}

\usage{
BICompare(phylo, t)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{t}{
the number of modalities to be tested}
 }

\value{
	BIC values for finding t modalities in the distance matrix of a tree and a control matrix, and a list of which nodes in the tree belong to each of t modalities.
}


\author{E Lewitus}


\examples{
data(Cetacea)
BICompare(Cetacea,5,method="gaussian")
}
