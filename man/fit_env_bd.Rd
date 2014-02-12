\name{fit_env_bd}
\alias{fit_env_bd}

\title{
Fit parameters for cladogenesis birth-death proces
}
\description{
This function helps to fit parameters of the birth-death model to a given phylogeny, by maximizing the likelihood function.
Outputs are the log-likelihood, the second order Akaike's Information Criterion, and the maximum likelihood estimates of the parameters of diversification. 
See notations in the PNAS 2011 paper.
}
\usage{
fit_env_bd(phylo, env_data, dof, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1,
       meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE,
       expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
       cond = "crown")
}

\arguments{
  \item{phylo}{
phylo is an object of type 'phylo' (see ape documentation)
}
  \item{env_data}{
Environmental data given as data frame with 2 columns
First column is time, second column is the envirionmental data (temperature for example)
Spline model is evaluated for the approximation of the envirionmental function
}
  \item{dof}{
Degree of freedom for the spline evalluation
}
  \item{tot_time}{
tot_time is the age of the phylogeny (crown age or stem age)
}
  \item{f.lamb}{
The \eqn{\lambda} function. This corresponds to the speciation rate function at any time.
The first argument is always time. The other arguments are the
  parameters of the variation, given as a single vector.
}
  \item{f.mu}{
The \eqn{\mu} function. This corresponds to the extinction rate function at any time.
The first argument is always time. The other arguments are the parameters of the variation, given as a single vector.
}

  \item{lamb_par}{
Initial values of lambda parameters
}
  \item{mu_par}{
Initial values of mu parameters
}
  \item{f}{
Probability that some extant species were sampled and not included
}
  \item{meth}{
The optimization method that should be used for the maximization of the likelihood function.
See \emph{optim} for more details.
}
  \item{cst.lamb}{
The flag fixes if the speciation rate is a constant function.
This parameter helps to perform numerical computations of the likelihood function.
}
  \item{cst.mu}{
The flag fixes if the extinction rate is a constant function.
This parameter helps to perform numerical computations of the likelihood function.
}
  \item{expo.lamb}{
The flag fixes if the speciation rate is an exponential function.
}
  \item{expo.mu}{
The flag fixes if the extinction rate is an exponential function.
}
  \item{fix.mu}{
If TRUE, then the mu parameter is fixed and will not be optimized
}
  \item{cond}{
cond is the conditioning to be used for the phylogeny.
Default value is cond=FALSE for no conditioning. Possibles values are:
\itemize{ \item FALSE for no conditioning;
\item "stem" is for conditioning on survival of the stem lineage (case when stem age is known);
\item "crown" is for conditioning on speciation and survival of 2 lineages (case when we work with crown ages);
}
}
}

\details{
Given a phylogeny, its associated likelihood function accounts at least one lineage surviving to present and being sampled.
Let \eqn{\Psi(s,t)} the probability that a lineage alive at time t leaves exactly one descendant
at time s (s < t).
Let \eqn{\Phi(s,t)} the probability that a lineage alive at time t has no descendant. The likelihood writes:
\deqn{L(t_{1},...,t_{n})=\frac{f^n\Psi(t_{2}, t_{1})\prod_{i=2}^{n} \Psi(s_{i,1}, t_{i}) \Psi(s_{i,2}, t_{i}) }{1-\Phi(t_{1})}}{see formula in Morlon et al. PNAS 2011}
where:
\itemize{\item f denotes the probability that some extant species were sampled and not included in the sample
\item \eqn{s_{i,1},s_{i,2}}{s[i,1],s[i,2]} denote the times at which the daughter lineages introduced at time \eqn{t_{i}}{t[i]} themselves branch.
}

This function fits this likelihood models. The function returns the the estimates
of the parameters of diversification, the log-likelihood and the second order Akaike's
Information Criterion, and the maximum likelihood.

The model likelihood is maximized using methods available in \code{\link[stats]{optim}}.
}
\value{
  \item{model}{The name of the used process}
  \item{LH}{The log-likelihood value}
  \item{aicc}{The second order Akaike's Information Criterion}
  \item{lamb_par}{Optimal lambda parameter}
  \item{mu_par}{Optimal mu parameters (if fix.mu is FALSE)}
}

\references{
Hélène Morlon, Matthew D. Potts, Joshua B. Plotkin: \emph{Inferring the Dynamics of Diversification: a Coalescent Approach}, PLoS Biology September 2010, Volume 8, Issue 9
}

\seealso{
\code{\link[RPANDA]{likelihood_bd}}
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
result <- fit_env_bd(Cetacea,env_data = InfTemp,dof=4, max(node.age(Cetacea)$ages),f.lamb,f.mu,
      lamb_par,mu_par,f=87/89, expo.lamb = TRUE, fix.mu=TRUE)
}
