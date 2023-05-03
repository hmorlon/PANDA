\name{div.models}
\alias{div.models}

\title{
  Model comparison of diversification models
}
\description{
  Applies a set of birth-death models to a phylogeny.
}
\usage{
  div.models(phylo, tot_time, f,
             backbone = F, spec_times = NULL, branch_times = NULL,
             models = c("BCST", "BCST_DCST", "BVAR",
                        "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR"),
             cond, verbose = T, n.max = NULL, rate.max = NULL)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
  }
  \item{f}{
numeric. The sampling fraction given as the number of species in the phylogeny over the number of species described in the taxonomy.
  }
  \item{backbone}{
  character. Allows to analyse a backbone. Default is FALSE and spec_times and branch_times are then ignored.
  Otherwise:
    \itemize{ 
      \item"stem.shift": for every shift, the probability of the speciation event at the stem age of the subclade is included in the likelihood of the backbone thanks to the argument spec_times.
      \item"crown.shift": for every shift, both the probability of the speciation event at the stem age of the subclade and the probability that the stem of the subclade survives to the crown age are included in the likelihood of the backbone thanks to the argument branch_times.
    }
  }
  \item{spec_times}{
  a numeric vector of the stem ages of subclades. Used only if backbone = "stem.shift". Default is NULL.
  }
  \item{branch_times}{
  a list of numeric vectors. Each vector contain the stem and crown ages of subclades (in this order). Used only if backbone = "crown.shift". Default is NULL.
  }
 \item{models}{
 a vector of character. Defines the set of birth-death models to applies e.g. BCST means pure-birth constant rate model, BCST_DVAR means birth constant rate and death variable rate model. Default is c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR") and applies all combinations of constant or variable rates for speciation and extinction. Time dependency is only exponential.
  }
  \item{cond}{
    conditioning to use to fit the model:
    \itemize{
      \item FALSE: no conditioning (not recommended);
      \item "stem": conditioning on the survival of the stem lineage (used when the stem age is known, in this case tot_time should be the stem age);
      \item "crown" (default): conditioning on a speciation event at the crown age and on survival of the two daugther lineages (used when the stem age is not known, in this case tot_time should be the crown age).
    }
  }
  \item{verbose}{
  bolean. Wether to print model names and AICc values during the calculation.
  }
  \item{rate.max}{
  numeric. Set a limit of diversificaton rates in terms of rate values.
  }
  \item{n.max}{
  numeric. Set a limit of diversificaton rates in terms of diversity estimates with the deterministic approach.
  }  
}
\details{
  Parameters of birth-death models are defined backward in time such as a positive alpha corresponds to a speciation rate decreasing through time from the past to the present.
}
\value{
  A data.frame with number of parameters, likelihood, AICc and parameter values for all models.
}
\references{
  Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (2023). Estimating clade-specific diversification rates and palaeodiversity dynamics from reconstructed phylogenies. BioRχiv.
}
\author{Nathan Mazet}

\examples{
data("Cetacea")
res <- div.models(Cetacea, tot_time = max(node.age(Cetacea)$ages),
                  f = 87/89, cond = "crown")
}