\name{plot_BICompare}
\alias{plot_BICompare}

\title{
Display modalities on a phylogeny.
}
\description{
Plot a phylogeny with branches colored according to modalities
}
\usage{
plot_BICompare(phylo,BICompare)
}

\arguments{
	 \item{phylo}{
an object of type 'phylo' (see ape documentation)}
  \item{BICompare}{
an object of class 'BICompare', output of the 'BICompare' function
}
	}


\value{
a plot of the phylogeny with branches colored according to which modalities they belong to.
}

\author{E Lewitus}

\references{
	Lewitus, E., Morlon, H., A graph-theoretical approach to characterising and comparing phylogenies (submitted) }


\seealso{
\code{\link{BICompare}}
}

\examples{

data(Cetacea)
result <- BICompare(Cetacea,5)
#plot_BICompare(Cetacea,result)
}
