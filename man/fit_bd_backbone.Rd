\name{fit_bd_backbone}
\alias{fit_bd_backbone}

\title{
Maximum likelihood fit of the general birth-death model (backbone)
}
\description{
Fits the birth-death model with potentially time-varying rates and potentially missing extant species to a phylogeny, by maximum likelihood. Notations follow Morlon et al. PNAS 2011. Modified version of fit_bd for backbones.
}
\usage{
fit_bd_backbone(phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1,
                backbone, spec_times, branch_times,
                meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE,
                expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
                dt=1e-3, cond = "crown", model)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
  \item{f.lamb}{
a function specifying the hypothesized functional form (e.g. constant, linear, exponential, etc.) of the variation of the speciation rate \eqn{\lambda} with time. Any functional form may be used.
This function has two arguments: the first argument is time; the second argument is a numeric vector of the parameters of the time-variation (to be estimated).
  }
  \item{f.mu}{
a function specifying the hypothesized functional form (e.g. constant, linear, exponential, etc.) of the variation of the extinction rate \eqn{\mu} with time. Any functional form may be used.
This function has two arguments: the first argument is time; the second argument is a numeric vector of the parameters of the time-variation (to be estimated).
  }
  \item{lamb_par}{
a numeric vector of initial values for the parameters of f.lamb to be estimated (these values are used by the optimization algorithm). The length of this vector is used to compute the total number of parameters in the model, so to fit a model with constant speciation rate (for example), lamb_par should be a vector of length 1. Otherwise aic values will be wrong.}
  \item{mu_par}{
a numeric vector of initial values for the parameters of f.mu to be estimated (these values are used by the optimization algorithm). The length of this vector is used to compute the total number of parameters in the model, so to fit a model without extinction (for example), mu_par should be empty (vector of length 0). Otherwise aic values will be wrong.
  }
  \item{f}{
the fraction of extant species included in the phylogeny
  }
  \item{backbone}{
    character. Allows to analyse a backbone. Default is FALSE and spec_times and branch_times are then ignored.
    Otherwise
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
  \item{meth}{
optimization to use to maximize the likelihood function, see \emph{optim} for more details.
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
  \item{fix.mu}{
 logical: if set to TRUE, the extinction rate \eqn{\mu} is fixed and will not be optimized.
  }
  \item{dt}{
the default value is 0. In this case, integrals in the likelihood are computed using R "integrate" function, which can be quite slow. If a positive dt is given as argument, integrals are computed using a piece-wise contant approximation, and dt represents the length of the intervals on which functions are assumed to be constant. For an exponential dependency of the speciation rate with time, we found that dt=1e-3 gives a good trade-off between precision and computation time.
  }
  \item{cond}{
    conditioning to use to fit the model:
    \itemize{ \item FALSE: no conditioning (not recommended);
      \item "stem": conditioning on the survival of the stem lineage (use when the stem age is known, in this case tot_time should be the stem age);
      \item "crown" (default): conditioning on a speciation event at the crown age and survival of the two daugther lineages (use when the stem age is not known, in this case tot_time should be the crown age).
    }
  }
  \item{model}{
  character. The model name as defined in the function div.models.
  }
}

\value{
	a list with the following components
  \item{model}{the name of the fitted model}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{lamb_par}{a numeric vector of estimated f.lamb parameters, in the same order as defined in f.lamb}
  \item{mu_par}{a numeric vector of estimated f.mu parameters, in the same order as defined in f.mu (if fix.mu is FALSE)}
}
\details{
The lengths of lamb_par and mu_par are used to compute the total number of parameters in the model, so to fit a model with constant speciation rate (for example), lamb_par should be a vector of length 1. Otherwise aic values will be wrong. In the f.lamb and f.mu functions, the first argument (time) runs from the present to the past. Hence, if the parameter controlling the variation of \eqn{\lambda} with time is estimated to be positive (for example), this means that the speciation rate decreases from past to present. Note that abs(f.lamb) and abs(f.mu) are used in the likelihood computation as speciation and extinction rates should always be positive. A consequence of 
this is that negative speciation/extinction rates estimates can be returned. They should be interpreted in absolute terms. See Morlon et al. 2020 for a more detailed explanation.}
\references{
	Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
	
	Morlon, H. (2014) Phylogenetic approaches for studying diversification, Eco Lett 17:508-525
	Morlon, H., Rolland, J. and Condamine, F. (2020) Response to Technical Comment ‘A cautionary note for users of linear diversification dependencies’, Eco Lett
	Mazet, N., Morlon, H., Fabre, P., Condamine, F.L., (2023). Estimating clade‐specific diversification rates and palaeodiversity dynamics from reconstructed phylogenies. \emph{Methods in Ecology and in Evolution} 14, 2575–2591. https://doi.org/10.1111/2041-210X.14195
}
\author{Hélène Morlon, Nathan Mazet}
\seealso{
\code{\link{plot_fit_bd}}, \code{\link{plot_dtt}}, \code{\link{likelihood_bd}, \link{fit_env}}
}
\examples{
# Some examples may take a little bit of time. Be patient!
data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
# Fit the pure birth model (no extinction) with a constant speciation rate
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par<-c(0.09)
mu_par<-c()
#result_cst <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,
#                     f=87/89,cst.lamb=TRUE,fix.mu=TRUE,dt=1e-3)
#result_cst$model <- "pure birth with constant speciation rate"
# Fit the pure birth model (no extinction) with exponential variation
# of the speciation rate with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
#result_exp <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,
#                     f=87/89,expo.lamb=TRUE,fix.mu=TRUE,dt=1e-3)
#result_exp$model <- "pure birth with exponential variation in speciation rate"
# Fit the pure birth model (no extinction) with linear variation of
# the speciation rate with time
f.lamb <-function(t,y){abs(y[1] + y[2] * t)}
# alternative formulation that can be used depending on the choice made to avoid negative rates: 
# f.lamb <-function(t,y){pmax(0,y[1] + y[2] * t)}, see Morlon et al. (2020)
f.mu<-function(t,y){0}
lamb_par<-c(0.09, 0.001)
mu_par<-c()
#result_lin <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89,fix.mu=TRUE,dt=1e-3)
#result_lin$model <- "pure birth with linear variation in speciation rate"
# Fit a birth-death model with exponential variation of the speciation
# rate with time and constant extinction
f.lamb<-function(t,y){y[1] * exp(y[2] * t)}
f.mu <-function(t,y){y[1]}
lamb_par <- c(0.05, 0.01)
mu_par <-c(0.005)
#result_bexp_dcst <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,
#                           f=87/89,expo.lamb=TRUE,cst.mu=TRUE,dt=1e-3)
#result_bexp_dcst$model <- "birth-death with exponential variation in speciation rate
#                           and constant extinction"
# Find the best model
#index <- which.min(c(result_cst$aicc, result_exp$aicc, result_lin$aicc,result_bexp_dcst$aicc))
#rbind(result_cst, result_exp, result_lin, result_bexp_dcst)[index,]
}