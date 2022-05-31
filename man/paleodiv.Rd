\name{paleodiv}
\alias{paleodiv}

\title{
  Paleodiversity through time
}
\description{
  Calculates paleodiversity through time from shift.estimates output with the deterministic approach.
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
    a data.frame containing a database of monophyletic groups for which potential shifts can be investigated. This database should be based on taxonomy, ecology or traits and contain a column named "Species" with species name as in phylo.
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

# testing family and parvorder levels
taxo_cetacea_no_genus <- taxo[names(taxo) != "Genus"]

# with f_df_cetacea the output from get.sampling.fractions()
# with shifts_cetacea the output from shift.estimates()
# use of paleodiv
paleodiversity <- paleodiv(phylo = Cetacea,
                           data = taxo_cetacea_no_genus,
                           sampling.fractions = f_df_cetacea,
                           shift.res = shifts_cetacea,
                           combi = 1, split.div = F)
}
}