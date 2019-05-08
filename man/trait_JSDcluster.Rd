\name{trait_JSDcluster}
\alias{trait_JSDcluster}

\title{
Clustering on the Jensen-Shannon distance between phylogenies
}

\description{
Computes the Jensen-Shannon distance metric between spectral density profiles of phylogenetic trait data and clusters on those distances. 
}

\usage{
trait_JSDcluster(phylo,mat,plot=F)
}

\arguments{
 \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{mat}{
  	a matrix of trait data with one trait per column and rows aligned to phylo tips}
  \item{plot}{
  	plot hierarchical cluster in a new window}	
   }


\value{
	plots a heatmap and hierarchical cluster with bootstrap support (>0.9) and outputs results of the k-medoids clustering on the optimal number of clusters in the form of a list with the following components
	\item{clusters}{a list with the following components: size, max_diss, av_diss, diameter, and separation}
	\item{J-S matrix}{a matrix providing the Jensen-Shannon distance values between pairs of phylogenetic trait data}
	\item{cluster assignment}{a table	that lists for each trait its cluster assignment and silhouete width}
	}
	

\author{E Lewitus}

\references{
	Lewitus, E., Morlon, H. (2015) Characterizing and comparing phylogenies from their Laplacian spectrum, Systematic Biology, 65(3):495-507}



\seealso{
\code{\link{trait_spectR}}
}


\examples{
tree<-data(Cetacea)
mat<-matrix(rexp(100, rate=.1), nrow=length(tree$tip.label),ncol=20)
trait_JSDcluster(tree,mat)
}
