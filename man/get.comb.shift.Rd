\name{get.comb.shift}
\alias{get.comb.shift}

\title{
  Combinations of shifts of diversification.
}
\description{
  Provides all the combinations of nodes of a phylogeny where shifts of diversification can be tested.
}
\usage{
  get.comb.shift(phylo, data, sampling.fractions, clade.size = 5, Ncores = 1)
}

\arguments{
  \item{phylo}{
    an object of type 'phylo' (see ape documentation)
  }
  \item{data}{
    an object of class "data.frame" containing a database of monophyletic groups for which potential shift can be investigated. This database should be based on taxonomy, ecology or traits and contain a column named "Species" with species name as in phylo.
  }
  \item{sampling.fractions}{
    the output resulting from get.sampling.fractions.
  }
  \item{clade.size}{
    numeric. Define the minimum number of species in a subgroup. Default is 5.
  }
  \item{Ncores}{
    numeric. Define the number of CPU to use for parallelizing the computation of combinations.
  }
}

\value{
  a list of combinations of subclades. The name of each element of list corresponds to the nodes of the subclades composing the combination separated by a point. If an combination cannot be splitted, its item in the list is NULL. Otherwise, it contains the internal nodes where the backbone can be splitted.
  
}
\details{
  More to come soon...

}
\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
	
}
\author{Nathan Mazet}
\seealso{
\code{\link{get.sampling.fractions}}, \code{\link{shift.estimates}}}
}
\examples{
# More to come soon...

}