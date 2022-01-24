\name{paleodiv}
\alias{div.rates}

\title{
  Daleodiversity through time
}
\description{
  Calculate paleodiversity through time from shift.estimates output.
}
\usage{
  paleodiv(phylo, data, sampling.fractions, shift.res,
           backbone.option = "crown.shift", combi = 1, split.div = F)
}

\arguments{
  \item{phylo}{
    an object of type 'phylo' (see ape documentation)
  }
  \item{data}{
    an object of class "data.frame" containing a database of monophyletic groups for which potential shift can be investigated. This database should be based on taxonomy, ecology or traits and contain a column named "Species" with species name as in phy.
  }
  \item{shift.res}{
    the output resulting from shift.estimates.
  }
  \item{backbone.option}{
    the output resulting from shift.estimates.
  }
  
  \item{combi}{
    numeric. The combination of shifts defined by its rank in the global comparison.
  }
  \item{split.div}{
    bolean. Specifies if paleodiversity should be return by parts
  }
}

\value{
  a matrix with as many rows as parts in the combination and as many column as million years from the root to the present.
  
}
\details{
  ...
  
}
\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
	
}
\author{Nathan Mazet}
\seealso{
\code{\link{shift.estimates}, \link{apply_prob_dtt}}
}
\examples{
# Examples are coming...

}