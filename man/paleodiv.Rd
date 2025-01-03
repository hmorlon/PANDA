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
           time.interval = 1, split.div = FALSE)
}

\arguments{
  \item{phylo}{
    an object of type 'phylo' (see ape documentation)
  }
  \item{data}{
    a data.frame containing a database of monophyletic groups for which potential shifts can be investigated. This database should be based on taxonomy, ecology or traits and contain a column named "Species" with species name as in phylo.
  }
  \item{sampling.fractions}{
    the output resulting from get.sampling.fractions.
  }
  \item{shift.res}{
    the output resulting from shift.estimates.
  }
  \item{backbone.option}{
    type of the backbone analysis:
    \itemize{ 
      \item"stem.shift": paleodiversity dynamics are calculated from the stem age for subclades.
      \item"crown.shift": paleodiversity dynamics are calculated from the crown age for subclades.
    }
  }
  
  \item{combi}{
    numeric. The combination of shifts defined by its rank in the global comparison.
  }
   \item{time.interval}{
  numeric. Define the time interval (in million years) at which paleodiversity values are calculated. Default is 1 for a value at each million year.
  }
  \item{split.div}{
    bolean. Specifies if paleodiversity should be plitted by parts of the selected combination (TRUE) or not. 
  }
}

\value{
  If split.div = FALSE, paleodiversity dynamics are returned in a matrix with as many rows as parts in the selected combination and as many column as million years from the root to the present. If spit.div = TRUE, global paleodiversity dynamic is returned as a vector with a value per million year.
}

\references{
  Mazet, N., Morlon, H., Fabre, P., Condamine, F.L., (2023). Estimating clade‐specific diversification rates and palaeodiversity dynamics from reconstructed phylogenies. \emph{Methods in Ecology and in Evolution} 14, 2575–2591. https://doi.org/10.1111/2041-210X.14195
}
\author{Nathan Mazet}
\seealso{
\code{\link{shift.estimates}, \link{apply_prob_dtt}}
}
\examples{
# loading data
data("Cetacea")
data("taxo_cetacea")
data("shifts_cetacea")

# no shifts tested at genus level
taxo_cetacea_no_genus <- taxo_cetacea[names(taxo_cetacea) != "Genus"]
f_cetacea <- get.sampling.fractions(phylo = Cetacea, lad = FALSE,
                                    data = taxo_cetacea_no_genus,
                                    plot = TRUE, cex = 0.3)
# use of paleodiv
paleodiversity <- paleodiv(phylo = Cetacea,
                           data = taxo_cetacea_no_genus,
                           sampling.fractions = f_cetacea,
                           shift.res = shifts_cetacea,
                           combi = 1, split.div = FALSE)
}