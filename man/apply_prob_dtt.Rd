\name{apply_prob_dtt}
\alias{apply_prob_dtt}

\title{
  Apply prob_dtt to output from shift.estimtes.
}
\description{
  Applies a set of birth-death models to a phylogeny.
}
\usage{
apply_prob_dtt(phylo, data, sampling.fractions, shift.res,
               combi = 1, backbone.option = "crown.shift")
}
\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
 \item{data}{
    an object of class "data.frame" containing a database of monophyletic groups for which potential shift can be investigated. This database should be based on taxonomy, ecology or traits and contain a column named "Species" with species name as in phy.
  }
  \item{sampling.fractions}{
    the output resulting from get.sampling.fractions.
  }
  \item{shift.res}{
    the output resulting from shift.estimates.
  }
  \item{backbone.option}{
    the output resulting from shift.estimates.
  }
  \item{combi}{
    numeric. The combination of shifts defined by its rank in the global comparison.
  }
}
\details{
  This funcion calls the function prob_dtt to calculate paleodiversity dynamic with the probabilistic 
}

\value{
  A list of results from prob_dtt for subclades and backbone(s).
}
\references{
Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record. Proc. Nat. Acad. Sci. 108: 16327-16332.

Billaud, O., Moen, D. S., Parsons, T. L., Morlon, H. (under review) Estimating Diversity Through Time using Molecular Phylogenies: Old and Species-Poor Frog Families are the Remnants of a Diverse Past.

Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep)
  
}

\author{Nathan Mazet}

\seealso{
\code{\link{fit_bd}, \link{plot_prob_dtt}, \link{prob_dtt}}
}

\examples{
# More to come soon...
}