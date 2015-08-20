\name{JSDtree}
\alias{JSDtree}

\title{
Jensen-Shannon distance between phylogenies
}

\description{
Computes the Jensen-Shannon distance metric between spectral density profiles of phylogenies. 
}

\usage{
JSDtree(phylo,meth=c("standard"))
}

\arguments{
  \item{phylo}{
a list of objects of type 'phylo' (see ape documentation)}
   \item{meth}{
the method used to compute the spectral density, which can either be "standard", "normal1", or "normal2". if set to "normal1", computes the spectral density normalized to the degree matrix. if set to "normal2", computes the spectral density normalized to the number of eigenvalues. if set to "standard", computes the unnormalized version of the spectral density (see the associated paper for an explanation)}
   }


\value{
	a matrix providing the Jensen-Shannon distance values between phylogeny pairs
		}

\author{E Lewitus}

\references{
Lewitus E., Morlon H. A graph-theoretical approach to characterising and comparing phylogenies (submitted)
}

\seealso{
\code{\link{JSDtree_cluster}}, \code{\link{spectR}}, \code{\link{BICompare}}
}


\examples{
trees<-TESS::tess.sim.age(n=20,age=10,0.15,0.05,MRCA=TRUE)
JSDtree(trees)
}
