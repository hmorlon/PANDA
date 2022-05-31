\name{apply_prob_dtt}
\alias{apply_prob_dtt}

\title{
  Calculates paleodiversity dynamics with the probabilistic approach.
}
\description{
  Applies prob_dtt to outputs from shift.estimates.
}
\usage{
apply_prob_dtt(phylo, data, sampling.fractions, shift.res,
               combi = 1, backbone.option = "crown.shift",
               m = NULL)
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
      \item"stem.shift": the stems of subclades are included in subclade analyses;
      \item"crown.shift": the stems of subclades are included in the backbone analysis (default).
    }
  }
  \item{combi}{
    numeric. The combination of shifts defined by its rank in the global comparison.
  }
  \item{m}{
    NULL or numeric. The set of maximum values for m ranges. Should be as long as the number of parts in the combinaison. Default is NULL (see details).
  }
  
}
\details{
  This funcion calls the function prob_dtt() to calculate paleodiversity dynamics with the probabilistic approach for the different parts of a combination of diversification shifts.
  
  As explained in Billaud et al. (2020), all the sum of probabilities per million year must be equal to 1. However, it can be difficult to reach 1 for groups showing a paleodiversity decline because the range of paleodiversity over which we need to calculate the probabilities can be very large. To circumvent this issue, apply_prob_dtt() set the range of the paleodiversity to the maximum of the deterministic estimate from the function `paleodiv()` and successively multiplies this maximum by 2, 3, 5, 7 and 10 until the sums of probabilities for each million year reach a minimum of 95\%. In few cases, this value of 95\% is not reached for few million years. In this case, it might come from an extremely high range of m and maximum values can be manually set up with the argument m.
}

\value{
  A list of results from prob_dtt() for subclades and backbone(s).
}
\references{
Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record. Proc. Nat. Acad. Sci. 108: 16327-16332.

Billaud, O., Moen, D.S., Parsons, T.L., Morlon, H., (2020). Estimating Diversity Through Time Using Molecular Phylogenies: Old and Species-Poor Frog Families are the Remnants of a Diverse Past. Systematic Biology 69, 363–383.

Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep)
}

\author{Nathan Mazet}

\seealso{
\code{\link{fit_bd}, \link{plot_prob_dtt}, \link{prob_dtt}}
}

\examples{
\dontrun{
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

# use of apply_prob_dtt()
rates <- apply_prob_dtt(phylo = Cetacea, shift.res = shifts_cetacea,
                        sampling.fractions = f_df_cetacea, combi = 1)
  
  }
}