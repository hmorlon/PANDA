\name{div.rates}
\alias{div.rates}

\title{
  Diversification rates through time
}
\description{
  Calculate diversification rates through time from shift.estimates output.
}
\usage{
  div.rates(phylo, shift.res, combi = 1,
  part = "backbone", backbone.option = "backbone2")
}

\arguments{
  \item{phylo}{
    an object of type 'phylo' (see ape documentation)
  }
  \item{shift.res}{
    the output resulting from shift.estimates.
  }
  \item{combi}{
    numeric. The combination of shifts defined by its rank in the global comparison.
  }
  \item{parts}{
    character. Specifies for which parts of the combination diversification rates has to be calculated. Default is "backbone" and provide only the backbone rate. Can be "all" for all the parts of a combination or "subclades" for subclades only.
  }
}

\value{
  a list of matrix with two rows (speciation and extinction) and as many column as million years from the root to the present.
  
}

\details{
  ...
}

\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
}

\author{Nathan Mazet}

\seealso{
\code{\link{shift.estimates}}
}
\examples{
# More coming soon...

}