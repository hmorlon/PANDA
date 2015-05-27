\name{BICompare}
\alias{BICompare}

\title{
Identify modalities in a phylogeny
}

\description{
Computes the BIC values for a specified number of modalities in the distance matrix of a phylogenetic tree and a control matrix and identifies these modalities using k-means clustering.
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
	a list with the following components:
	\item{BIC_test}{
	BIC values for finding t modalities in the distance matrix of a tree and a control matrix}
	\item{clusters}{
	 a vector specifying which nodes in the tree belong to each of t modalities}
	 The function also returns a plot of the phylogeny with branches colored according to which modality it belongs to.
	 
}


\author{E Lewitus}


\examples{
data(Cetacea)
BICompare(Cetacea,5)
}
