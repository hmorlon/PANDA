\name{simul.comb.shift}
\alias{simul.comb.shift}

\title{
  Simulating trees from shift.estimates() results to test model adequacy
}
\description{
  Simulates trees with combination of shifts from shifts.estimates() output.
}
\usage{
  simul.comb.shift(n = 10000, phylo, sampling.fractions,
                    shift.res, combi = 1, clade.size = 5)
}
\arguments{
  \item{n}{
    numeric. Defines the number of simulations to generate (see Details).
  }
  \item{phylo}{
    an object of type 'phylo' (see ape documentation).
  }
  \item{sampling.fractions}{
    the output resulting from get.sampling.fractions.
  }
  \item{shift.res}{
    the output resulting from shift.estimates.
  }
  \item{combi}{
    numeric. Corresponds to the rank of the combination in the global comparison (shift.res$total).
  }
  \item{clade.size}{
    numeric. Defines the minimum number of species in a subgroup. Default is 5.
  }
}

\value{
  a list of simulated phylogenies as object of type phylo. Tips of subclades are named with the letters a, b, c, etc. while tips of backbones are named with letters z, y, etc. The empirical groups are sorted from the more recent to the older one (i.e. group a will be the more recent empirical subclade, etc.).
}

\details{
  Some combinations of shifts might be complex cases to simulate because the backbone needs to be rich enough to graft subclades. Some simulations will not satisfy this condition and will then be discarded. In consequence, the number of simulated phylogenies in the output will not be equal to n for complex simulations. This is why the value of n is high by default (n = 10000), to ensure to have enough simulations (around 500) to test the presence.

clade.size argument should be the same value for the whole procedure in the empirical case (same than for get.sampling.fraction and get.comb.shift).

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
data("taxo_cetacea")
data("shifts_cetacea")

# with the results from shifts.estimates()

all_posteriors_cetacea <- simul_posteriors(phylo = Cetacea,
                                           sampling.fractions = f_cetacea,
                                           shift.res = shifts_cetacea)
  }
}
