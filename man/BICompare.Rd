\name{BICompare}
\alias{BICompare}

\title{
Identify modalities in a phylogeny
}

\description{
Computes the BIC values for a specified number of modalities in the distance matrix of a phylogenetic tree and that of a random tree; identifies these modalities using k-means clustering.
}

\usage{
BICompare(phylo,t)
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
	BIC values for finding t modalities in the distance matrix of a tree and that of a random ("control") tree}
	\item{clusters}{
	 a vector specifying which nodes in the tree belong to each of t modalities}
	 \item{BSS/TSS}{
	 the ratio of between-cluster sum of squares over total sum of squares}
	 
}

\author{E Lewitus}

\references{
	Lewitus, E., Morlon, H., A graph-theoretical approach to characterising and comparing phylogenies (submitted) }


\seealso{
\code{\link{plot_BICompare}}, \code{\link{spectR}}, \code{\link{JSDtree}}
}


\examples{
data(Cetacea)
BICompare(Cetacea,5)
}
