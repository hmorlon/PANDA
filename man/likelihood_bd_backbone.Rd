\name{likelihood_bd_backbone}
\alias{likelihood_bd_backbone}

\title{
Likelihood of a phylogeny under the general birth-death model (backbone)
}
\description{
Computes the likelihood of a phylogeny under a birth-death model with potentially time-varying rates and potentially missing extant species. Notations follow Morlon et al. PNAS 2011. Modified version of likelihood_bd for backbones.
}
\usage{
likelihood_bd_backbone(phylo, tot_time, f, f.lamb, f.mu, 
                       backbone, spec_times, branch_times,
                       cst.lamb = FALSE, cst.mu = FALSE,
                       expo.lamb = FALSE, expo.mu = FALSE, dt=0, cond = "crown")
}

\arguments{
  \item{phylo}{  
an object of type 'phylo' (see ape documentation)
}
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
 \item{f.lamb}{
  a function specifying the time-variation of the speciation rate \eqn{\lambda}. This function as a single argument (time). Any function may be used.
}
  \item{f.mu}{
  a function specifying the time-variation of the speciation rate \eqn{\mu}. This function as a single argument (time). Any function may be used.
}
 \item{f}{
the fraction of extant species included in the phylogeny
}
  \item{backbone}{
  character. Allows to analyse a backbone. Default is NULL and spec_times and branch_times are then ignored.
  
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
  a list of numeric vectors. Each vector contains the stem and crown ages of subclades (in this order). Used only if backbone = "crown.shift". Default is NULL.
  }
 \item{cst.lamb}{
logical: should be set to TRUE only if f.lamb is constant (i.e. does not depend on time) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{cst.mu}{
logical: should be set to TRUE only if f.mu is constant (i.e. does not depend on time) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.lamb}{
logical: should be set to TRUE only if f.lamb is exponential to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.mu}{
logical: should be set to TRUE only if f.mu is exponential to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{dt}{
 the default value is 0. In this case, integrals in the likelihood are computed using R "integrate" function, which can be quite slow. If a positive dt is given as argument, integrals are computed using a piece-wise contant approximation, and dt represents the length of the intervals on which functions are assumed to be constant. For an exponential dependency of the speciation rate with time, we found that dt=1e-3 gives a good trade-off between precision and computation time.
}
  \item{cond}{
conditioning to use to fit the model:
\itemize{ \item FALSE: no conditioning (not recommended);
\item "stem": conditioning on the survival of the stem lineage (use when the stem age is known, in this case tot_time should be the stem age);
\item "crown" (default): conditioning on a speciation event at the crown age and survival of the 2 daugther lineages (use when the stem age is not known, in this case tot_time should be the crown age).
}
}
}
\details{
	When specifying f.lamb and f.mu, time runs from the present to the past (hence if the speciation rate decreases with time, f.lamb must be a positive function of time).
}

\value{
  the loglikelihood value of the phylogeny, given f.lamb and f.mu
}
\references{
	Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
	Mazet, N., Morlon, H., Fabre, P., Condamine, F.L., (2023). Estimating clade‐specific diversification rates and palaeodiversity dynamics from reconstructed phylogenies. \emph{Methods in Ecology and in Evolution} 14, 2575–2591. https://doi.org/10.1111/2041-210X.14195
}

\author{Hélène Morlon, Nathan Mazet}

\examples{
data(Cetacea)
tot_time <- max(node.age(Cetacea)$ages)
# Compute the likelihood for a pure birth model (no extinction) with
# an exponential variation of speciation rate with time
lamb_par <- c(0.1, 0.01)
f.lamb <- function(t){lamb_par[1] * exp(lamb_par[2] * t)}
f.mu <- function(t){0}
f <- 87/89
# same as likelihood_bd in this case
lh <- likelihood_bd_backbone(Cetacea, tot_time, f, f.lamb, f.mu, 
                             backbone = FALSE, spec_times = NULL, branch_times = NULL,
                             cst.mu = TRUE, expo.lamb = TRUE, dt = 1e-3)
}