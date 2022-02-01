\name{shift.estimates}
\alias{shift.estimates}

\title{
  Estimating shifts of diversification
}
\description{
  Applies models of diversification to each part of all combinations of shifts to detect the best combination of subclades and backbone(s).
}
\usage{
  shift.estimates(phylo, data, sampling.fractions, comb.shift,
                  models = c("BCST", "BCST_DCST", "BVAR",
                  "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR"),
                  backbone.option = "crown.shift", multi.backbone = F,
                  np.sub = 4, rate.max = NULL, n.max = NULL, Ncores = 1)
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
  \item{comb.shift}{
    the output resulting from get.comb.shift.
  }
  \item{models}{
    a vector of character that specifies the set of model of diversification to apply. Default is c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR").
  }
  \item{backbone.option}{
    type of the backbone analysis:
    \itemize{ 
      \item"stem.shift": the stems of subclades are included in subclade analyses;
      \item"crown.shift": the stems of subclades are included in the backbone analysis (Default).
    }
  }
  
  \item{multi.backbone}{
    can be either FALSE (default), TRUE or "all": 
    \itemize{
    \item{FALSE: only combinations with simple backbone will be analyzed.}
    \item{TRUE: only multiple backbones will be analyzed.}
    \item{"all", all combinations are analyzed.}
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
  \item{backbones}{a list with the model comparison for all subclades (see details)}
  \item{total}{the global comparison of combinations based on AICc}
  
}
\details{
  The output for backbone is a list of with three levels. The first level corresponds to the combinations of subclades. The second level corresponds to the potential splitting within the backbone if analysis with multiple backbones has been led. The third level is the list of dataframe summaryzing the model comparison for all backbone of this combination of subclade (see example section).

We recommand to remove "BVAR_DVAR" model from the models set and to lead the first analysis with multi.backbone = F.

clade.size argument should be the same value for the whole procedure (same than for get.sampling.fraction and get.comb.shift).

}
\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
	
}
\author{Nathan Mazet}
\seealso{
\code{\link{get.sampling.fractions}}, \code{\link{shift.estimates}}, code{\link{paleodiv}}
}
\examples{
# loading data

data("Cetacea")
data("taxo_cetacea")

# main procedure
taxo_cetacea_no_genus <- taxo[names(taxo) != "Genus"]
f_df_cetacea <- get.sampling.fractions(phylo = Cetacea, lad = F,
                                             data = taxo_cetacea_no_genus,
                                             plot = T, cex = 0.3)

comb.shift_cetacea <- get.comb.shift(phylo = Cetacea, data = taxo[names(taxo) != "Genus"],
                             sampling.fractions = f_df_cetacea, Ncores = 4)

shifts_cetacea <- shift.estimates(phylo = Cetacea, data = taxo_cetacea,
                          sampling.fractions = f_df_cetacea,
                          comb.shift = comb.shift_cetacea,
                          models = c("BCST","BCST_DCST","BVAR",
                          "BVAR_DCST","BCST_DVAR"),
                          backbone.option = "crown.shift", Ncores = 4)

}