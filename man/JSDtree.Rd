\name{JSDtree}
\alias{JSDtree}

\title{
The Jensen-Shannon distance metric between phylogenies
}

\description{
Computes the Jensen-Shannon distance metric between spectral density profiles of phylogenies.
}

\usage{
JSDtree(phylo,meth=c("standard"),alpha)
}

\arguments{
  \item{phylo}{
a list of objects of type 'phylo' (see ape documentation)}
   \item{meth}{
the method used to compute the spectral density, which can either be "standard", "normal1", or "normal2". if set to "normal1", compute the spectral density normalized to the degree matrix. if set to "normal2", compute the spectral density normalized to the number of eigenvalues. if set to "standard", compute the unnormalized version of the spectral density (see the associated paper for an explanation)}
   \item{alpha}{
   	the confidence value for demarcating clusters in the hierarchical clustering plot
   }
}


\value{
	a list with the following components:
	\item{JSD}{a matrix providing the Jensen-Shannon distance values between phylogeny pairs}
	\item{clusters}{the optimal number of clusters around medoids (see pamk documentation)}
	\item{cluster_assignments}{assignments of trees to clusters} 
	}


\references{
Lewitus E., Morlon H. (submitted) A graph-theoretical approach to characterising and comparing phylogenies
}

\author{E Lewitus}

\seealso{
\code{\link{plot_JSDtree}}, \code{\link{spectR}}, \code{\link{BICompare}}
}


\examples{
trees<-TESS::sim.globalBiDe.age(n=20,age=10,0.15,0.05,MRCA=TRUE)
JSDtree(trees)
}
