\name{get.comb.shift}
\alias{get.comb.shift}

\title{
  Combinations of shifts of diversification.
}
\description{
  Provides all the combinations of nodes of a phylogeny where shifts of diversification can be tested.
}
\usage{
  get.comb.shift(phylo, data, sampling.fractions,
                 clade.size = 5, Ncores = 1)
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
    numeric. Define the number of CPU cores to use for parallelizing the computation of combinations.
  }
}

\value{
  a vector of character summaryzing the combination of shifts as a concatenation of node IDs separated by "." or "/". Node IDs at the left of "/" correspond to shifts at the origin of subclades (monophyletic and ultrametric subtrees) while node IDs at the right of "/" correspond to shifts at the origin of backbone(s) (pruned trees).
  
}
\details{
  clade.size argument should be the same value for the whole procedure (same than for get.sampling.fraction and shift.estimates).
}
\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
	
}
\author{Nathan Mazet}
\seealso{
\code{\link{get.sampling.fractions}}, \code{\link{shift.estimates}}}

\examples{
  \dontrun{
  # loading data
  
  data("Cetacea")
  data("taxo_cetacea")
  
  # main procedure
  taxo_cetacea_no_genus <- taxo[names(taxo) != "Genus"]
  f_cetacea <- get.sampling.fractions(phylo = Cetacea, lad = F,
                                      data = taxo_cetacea_no_genus,
                                      plot = T, cex = 0.3)
  
  comb.shift_cetacea <- get.comb.shift(phylo = Cetacea,
                                       data = taxo_cetacea_no_genus,
                                       sampling.fractions = f_df_cetacea, Ncores = 4)
  }
}