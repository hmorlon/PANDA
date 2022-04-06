\name{div.rates}
\alias{div.rates}

\title{
  Diversification rates through time
}
\description{
  Calculates diversification rates through time from shift.estimates() output.
}
\usage{
  div.rates(phylo, shift.res, combi = 1,
            part = "backbone", backbone.option = "crown.shift")
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
  \item{backbone.option}{
    type of the backbone analysis:
    \itemize{ 
      \item"stem.shift": the stems of subclades are included in subclade analyses;
      \item"crown.shift": the stems of subclades are included in the backbone analysis (Default).
    }
  }
}

\value{
  a list of matrix with two rows (speciation and extinction) and as many columnS as million years from the root to the present.
}

\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
}

\author{Nathan Mazet}

\seealso{
\code{\link{shift.estimates}}
}
\examples{
\dontrun{
# loading data
data("Cetacea")

# with shifts_cetacea the output from shift.estimates()
rates <- div.rates(phylo = Cetacea, shift.res = shifts_cetacea,
                   combi = 1, part = "all")
}
}