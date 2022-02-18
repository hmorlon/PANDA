\name{paleodiv}
\alias{paleodiv}

\title{
  Paleodiversity through time
}
\description{
  Calculates paleodiversity through time from shift.estimates output with the determinsitic approach.
}
\usage{
  paleodiv(phylo, data, sampling.fractions, shift.res,
           backbone.option = "crown.shift", combi = 1,
           split.div = F)
}

\arguments{
  \item{phylo}{
    an object of type 'phylo' (see ape documentation)
  }
  \item{data}{
    a data.frame containing a database of monophyletic groups for which potential shift can be investigated. This database should be based on taxonomy, ecology or traits and contain a column named "Species" with species name as in phy.
  }
  \item{shift.res}{
    the output resulting from shift.estimates.
  }
  \item{backbone.option}{
    type of the backbone analysis:
    \itemize{ 
      \item"stem.shift": the stems of subclades are included in subclade analyses;
      \item"crown.shift": the stems of subclades are included in the backbone analysis (Default).
    }
  }
  
  \item{combi}{
    numeric. The combination of shifts defined by its rank in the global comparison.
  }
  \item{split.div}{
    bolean. Specifies if paleodiversity should be plitted by parts of the selected combination (TRUE) or not. 
  }
}

\value{
  If split.div = FALSE, paleodiversity dynamics are returned in a matrix with as many rows as parts in the selected combination and as many column as million years from the root to the present. If spit.div = TRUE, global paleodiversity dynamic is returned as a vector with a value per million year.
}

\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
	
}
\author{Nathan Mazet}
\seealso{
\code{\link{shift.estimates}, \link{apply_prob_dtt}}
}
\examples{
\dontrun{
# loading data
data("Cetacea")
data("taxo_cetacea")

# main procedure
taxo_cetacea_no_genus <- taxo[names(taxo) != "Genus"]
f_df_cetacea <- get.sampling.fractions(phylo = Cetacea, lad = F,
                                       data = taxo_cetacea_no_genus,
                                       plot = T, cex = 0.3)

comb.shift_cetacea <- get.comb.shift(phylo = Cetacea,
                                     data = taxo_cetacea_no_genus,
                                     sampling.fractions = f_df_cetacea,
                                     Ncores = 4)

shifts_cetacea <- shift.estimates(phylo = Cetacea,
                                  data = taxo_cetacea_no_genus,
                                  sampling.fractions = f_df_cetacea,
                                  comb.shift = comb.shift_cetacea,
                                  models = c("BCST","BCST_DCST","BVAR",
                                  "BVAR_DCST","BCST_DVAR"),
                                  Ncores = 4)

# use of paleodiv
paleodiversity <- paleodiv(phylo = Cetacea,
                           data = taxo_cetacea_no_genus,
                           sampling.fractions = f_df_cetacea,
                           shift.res = shifts,
                           combi = 1, split.div = F)
}
}