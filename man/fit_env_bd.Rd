\name{fit_env_bd}
\alias{fit_env_bd}

\title{
Maximum likelihood fit of the environmental birth-death model
}
\description{
Fits the environmental birth-death model with potentially missing extant species to a phylogeny, by maximum likelihood. Notations follow Morlon et al. PNAS 2011 and Condamine et al. ELE 2013.
}
\usage{
fit_env_bd(phylo, env_data, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1,
       meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE,
       expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
       cond = "crown",df)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{env_data}{
environmental data, given as a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).
}
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
  \item{f.lamb}{
a function specifying the hypothesized functional form of the variation of the speciation rate \eqn{\lambda} with time and the environmental variable. Any functional form may be used.
This function has three arguments: the first argument is time; the second argument is the environmental variable; the third arguement is a numeric vector of the parameters controlling the time and environmental variation (to be estimated).
}
  \item{f.mu}{
a function specifying the hypothesized functional form of the variation of the extinction rate \eqn{\mu} with time and the environmental variable. Any functional form may be used.
This function has three arguments: the first argument is time; the second argument is the environmental variable; the second argument is a numeric vector of the parameters controlling the time and environmental variation (to be estimated).
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
logical: should be set to TRUE only if f.lamb is constant (i.e. does not depend on time or the environmental variable) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{cst.mu}{
logical: should be set to TRUE only if f.mu is constant (i.e. does not depend on time or the environmental variable) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.lamb}{
logical: should be set to TRUE only if f.lamb is an exponential function of time (and does not depend on the environmental variable) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.mu}{
logical: should be set to TRUE only if f.mu is an exponential function of time (and does not depend on the environmental variable) to use analytical instead of numerical computation in order to reduce computation time.
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
\item{df}{
 the degree of freedom to use to define the spline. As a default, smooth.spline(env_data[,1], env_data[,2])$df may be used. See sm.spline for details.  
}

}


\value{
	a list with the following components
  \item{model}{the name of the used process}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{lamb_par}{a numeric vector of estimated f.lamb parameters, in the same order as defined in f.lamb}
  \item{mu_par}{a numeric vector of estimated f.mu parameters, in the same order as defined in f.mu (if fix.mu is FALSE)}
}

\details{
In the f.lamb and f.mu functions, time runs from the present to the past.
%The rate models here are function of environmental data fixed by
%env_data.
%Hence, if the parameter controlling the variation of \eqn{\lambda} with time is estimated to be positive (for example), this means that the speciation rate decreases from past to present.
}

\note{
	The speed of convergence of the fit might depend on the degree of freedom chosen to define the spline. 
}
\references{
Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
	
Condamine, F.L., Rolland, J., and Morlon, H. (2013) Macroevolutionary perspectives to environmental change, Eco Lett 16: 72-85

Morlon, H. (2014) Phylogenetic approaches for studying diversification, Eco Lett
}

\seealso{
\code{\link[PANDA]{fit_bd},\link[PANDA]{likelihood_bd}}
}

\examples{
data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
data(InfTemp)
df<-smooth.spline(InfTemp[,1], InfTemp[,2])$df

# fits a model with lambda varying as an exponential function of temperature and mu fixed to 0 (no extinction)

f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10, 0.01)
mu_par<-c()
result_exp <- fit_env_bd(Cetacea,InfTemp,tot_time,f.lamb,f.mu,
      lamb_par,mu_par,f=87/89, fix.mu=TRUE)
}
