\name{fit_env_bd}
\alias{fit_env_bd}

\title{
Maximum likelihood fit of the general birth-death model
}
\description{
Fits the environmental birth-death model with time-varying rates and potentially missing extant species to a phylogeny, by maximum likelihood. Notations follow Morlon et al. PNAS 2011.
}
\usage{
fit_env_bd(phylo, env_data, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1,
       meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE,
       expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
       cond = "crown")
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{env_data}{
Environmental data given as data frame with 2 columns
First column is time, second column is the environmental data (temperature for example)
Spline model is evaluated for the approximation of the environmental function
}
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
  \item{f.lamb}{
a function specifying the hypothesized functional form (e.g. constant, linear, exponential, etc.) of the variation of the speciation rate \eqn{\lambda} with time. Any functional form may be used.
This function as two arguments: the first argument is time; the second argument is a numeric vector of the parameters of the time-variation (to be estimated).
}
  \item{f.mu}{
a function specifying the hypothesized functional form (e.g. constant, linear, exponential, etc.) of the variation of the extinction rate \eqn{\mu} with time. Any functional form may be used.
This function as two arguments: the first argument is time; the second argument is a numeric vector of the parameters of the time-variation (to be estimated).
}

  \item{lamb_par}{
a numeric vector of initial values for the parameters of f.lamb to be estimated (these values are used by the optimization algorithm) 
}
  \item{mu_par}{
a numeric vector of initial values for the parameters of f.mu to be estimated (these values are used by the optimization algorithm) 
}
  \item{f}{
the fraction of extant species included in the phylogeny
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
logical: should be set to TRUE only if f.lamb is exponential (i.e. does not depend on time) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.mu}{
logical: should be set to TRUE only if f.mu is exponential (i.e. does not depend on time) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{fix.mu}{
 logical: if set to TRUE, the extinction rate \eqn{\mu} is fixed and will not be optimized.
}
  \item{cond}{
conditioning to use to fit the model:
\itemize{ \item FALSE: no conditioning (not recommended);
\item "stem": conditioning on the survival of the stem lineage (use when the stem age is known, in this case tot_time should be the stem age);
\item "crown" (default): conditioning on a speciation event at the crown age and survival of the 2 daugther lineages (use when the stem age is not known, in this case tot_time should be the crown age).
}
}
}


\value{
	a list with the following components
  \item{model}{The name of the used process}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{lamb_par}{a numeric vector of estimated f.lamb parameters, in the same order as defined in f.lamb}
  \item{mu_par}{a numeric vector of estimated f.mu parameters, in the same order as defined in f.mu (if fix.mu is FALSE)}
}

\details{
In the f.lamb and f.mu functions, the first argument (time) runs from the present to the past. 
The rate models here are function of environmental data fixed by env_data.
%Hence, if the parameter controlling the variation of \eqn{\lambda} with time is estimated to be positive (for example), this means that the speciation rate decreases from past to present.
}

\references{
	Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
}

\seealso{
\code{\link[PANDA]{likelihood_bd}}
}

\examples{
# Phylogenetic data
data(Cetacea)
# Environmental data: temperature
data(InfTemp)
#first model
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
# initial points for optimization
lamb_par<-c(0.05, 0.01)
mu_par<-c()
# fit the model
result <- fit_env_bd(Cetacea,env_data = InfTemp, max(node.age(Cetacea)$ages),f.lamb,f.mu,
      lamb_par,mu_par,f=87/89, expo.lamb = TRUE, fix.mu=TRUE)
}
