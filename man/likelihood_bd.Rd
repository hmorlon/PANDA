\name{likelihood_bd}
\alias{likelihood_bd}

\title{
Compute the likelihood of a phylogeny under the general birth-death model
}
\description{
Computes the likelihood of a phylogeny under a birth-death model with potentially time-varying rates and potentially missing extant species. Notations follow Morlon et al. PNAS 2011.
}
\usage{
likelihood_bd(phylo, tot_time, f.lamb, f.mu, f, cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE,
expo.mu = FALSE, cond = "crown")
}

\arguments{
  \item{phylo}{  
an object of type 'phylo' (see ape documentation)
}
  
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}

 \item{f.lamb}{
a function specifying the time-variation of the speciation rate \eqn{\lambda}. This function as a single argument (time).
}
  \item{f.mu}{
  	a function specifying the time-variation of the speciation rate \eqn{\mu}. This function as a single argument (time).
}

 \item{f}{
the fraction of extant species included in the phylogeny
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
%Given a phylogeny, its associated likelihood function accounts at least one lineage surviving to present and being sampled.
%Let \eqn{\Psi(s,t)} the probability that a lineage alive at time t leaves exactly one descendant
%at time s (s < t).
%Let \eqn{\Phi(s,t)} the probability that a lineage alive at time t has no descendant. Likelihood becomes:
%\eqn{L(t_{1},...,t_{n})=\frac{f^n\Psi(t_{2}, t_{1})\prod_{i=2}^{n} \Psi(s_{i,1}, t_{i}) \Psi(s_{i,2}, t_{i}) }{1-\Phi(t_{1})}}{L(t[1],...,t[n])=f^n(\Psi(t[2], t[1])/(1-\Phi(t[1])))*\prod_{i=2}^{n} \Psi(s[i,1], t[i]) \Psi(s[i,2], t[i])}
%where:
%%\itemize{
%\item f denotes the probability that some extant species were sampled and not included in the sample
%\item \eqn{s_{i,1},s_{i,2}}{s[i,1],s[i,2]} denote the times at which the daughter lineages introduced at time \eqn{t_{i}}{t[i]} themselves branch.}
}

\value{
  \item{loglikelihood}{the loglikelihood value of the phylogeny given f.lamb and f.mu}
}
\references{
	Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
}

\examples{
data(Cetacea)
tot_time <- max(node.age(Cetacea)$ages)
#compute the likelihood for a pure birth model (no extinction) with an exponential variation of speciation rate with time
lamb_par <- c(0.1, 0.01)
%# lambda(t) and mu(t)
f.lamb <- function(x){lamb_par[1] * exp(lamb_par[2] * x)}
f.mu <- function(x){0}
f <- 87/89
%# model with none extinction & exponential speciation variation
likelihood_bd(Cetacea,tot_time,f.lamb,f.mu,f,cst.mu=TRUE,expo.lamb=TRUE)
}
