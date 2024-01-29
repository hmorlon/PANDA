\name{shift.estimates}
\alias{shift.estimates}

\title{
  Estimating clade-shifts of diversification
}
\description{
  Applies models of diversification to each part of all combinations of shifts to detect the best combination of subclades and backbone(s).
}
\usage{
  shift.estimates(phylo, data, sampling.fractions, comb.shift,
                  models = c("BCST", "BCST_DCST", "BVAR",
                  "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR"),
                  backbone.option = "crown.shift",
                  multi.backbone = F, np.sub = 4,
                  rate.max = NULL, n.max = NULL, Ncores = 1)
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
  \item{comb.shift}{
    the output resulting from get.comb.shift.
  }
  \item{models}{
    a vector of character that specifies the set of models of diversification to apply. Default is c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR").
  }
  \item{backbone.option}{
    type of the backbone analysis:
    \itemize{ 
      \item"stem.shift": for every shift, the probability of the speciation event at the stem age of the subclade is included in the likelihood of the backbone thanks to the argument spec_times.
      \item"crown.shift": for every shift, both the probability of the speciation event at the stem age of the subclade and the probability that the stem of the subclade survives to the crown age are included in the likelihood of the backbone thanks to the argument branch_times.
    }
  }
  \item{multi.backbone}{
    can be either FALSE (default), TRUE or "all": 
    \itemize{
    \item{FALSE: only combinations with simple backbone will be analyzed.}
    \item{TRUE: only combinations with multiple backbones will be analyzed.}
    \item{"all": all combinations are analyzed.}
    }
  }
  \item{np.sub}{
    Defines the set of models to apply to subclade based on the number of parameters. By default np.sub = 4 and all models from argument models will be applied. If np.sub = 3, the more complex model "BVAR_DVAR" is excluded. If np.sub = 2, the set of models is reduced to "BCST", "BCST_DCST" and "BVAR" models. np.sub = "no_extinction" only applies "BCST" and "BVAR" models.
  }
  \item{rate.max}{
    numeric. Define a maximum value for diversification rate through time.
  }
  \item{n.max}{
    numeric. Define a maximum value for diversity through time.
  }
  \item{Ncores}{
    numeric. Define the number of CPU cores to use for parallelizing the computation of combinations.
  }
}

\value{
  a list with the following components
  \item{whole_tree}{a data.frame with the model comparison for the whole tree}
  \item{subclades}{a list of dataframes summaryzing the model comparison for all subclades (same format than div.models outputs)}
  \item{backbones}{a list with the model comparison for all backbones (see details)}
  \item{total}{the global comparison of combinations based on AICc}
  
}
\details{
  The output for backbone is a list in which each element corresponds to the backbone model comparisons of a combination. This element contains a list with one table of model comparison per backbone.

We recommand to remove "BVAR_DVAR" model from the models set and to lead the first analysis with multi.backbone = F to limit the number of combination.

clade.size argument should be the same value for the whole procedure (same than for get.sampling.fraction and get.comb.shift).

}
\references{
	Mazet, N., Morlon, H., Fabre, P., Condamine, F.L., (2023). Estimating clade‐specific diversification rates and palaeodiversity dynamics from reconstructed phylogenies. \emph{Methods in Ecology and in Evolution} 14, 2575–2591. https://doi.org/10.1111/2041-210X.14195
	
}
\author{Nathan Mazet}
\seealso{
\code{\link{get.sampling.fractions}}, \code{\link{shift.estimates}}, \code{\link{paleodiv}}
}
\examples{
\donttest{
# loading data
data("Cetacea")
data("taxo_cetacea")

# whole procedure
taxo_cetacea_no_genus <- taxo_cetacea[names(taxo_cetacea) != "Genus"]
f_cetacea <- get.sampling.fractions(phylo = Cetacea, lad = FALSE,
                                    data = taxo_cetacea_no_genus,
                                    plot = TRUE, cex = 0.3)

comb.shift_cetacea <- get.comb.shift(phylo = Cetacea,
                                     data = taxo_cetacea_no_genus,
                                     sampling.fractions = f_cetacea,
                                     Ncores = 4)
                                     
shifts_cetacea <- shift.estimates(phylo = Cetacea,
                                  data = taxo_cetacea_no_genus,
                                  sampling.fractions = f_cetacea,
                                  comb.shift = comb.shift_cetacea,
                                  models = c("BCST","BCST_DCST","BVAR",
                                             "BVAR_DCST","BCST_DVAR"),
                                  backbone.option = "crown.shift",
                                  Ncores = 4)
  }
}
