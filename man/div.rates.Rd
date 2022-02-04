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
# loading data

data("Cetacea")
data("taxo_cetacea")

# main procedure
f_df_cetacea <- get.sampling.fractions(phylo = Cetacea, lad = F,
                                             data = taxo[names(taxo) != "Genus"],
                                             plot = T, cex = 0.3)

comb.shift_cetacea <- get.comb.shift(phylo = Cetacea, data = taxo[names(taxo) != "Genus"],
                             sampling.fractions = f_df_cetacea, Ncores = 4)

shifts_cetacea <- shift.estimates(phylo = Cetacea, data = taxo_cetacea,
                          sampling.fractions = f_df_cetacea,
                          comb.shift = comb.shift_cetacea,
                          models = c("BCST","BCST_DCST","BVAR",
                          "BVAR_DCST","BCST_DVAR"),
                          backbone.option = "crown.shift", Ncores = 4)

# use of div.rates
rates <- div.rates(phylo = Cetacea, shift.res = shifts, combi = 1, part = "all")

}