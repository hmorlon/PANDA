\name{likelihood_bd}
\alias{likelihood_bd}

\title{
likelihood of a phylogeny conditioned to the presence of at least one descendant in the sample
}
\description{
The function computes the likelihood of a phylogeny of n species using a birth-death model of cladogenesis.
The cladogenesis process assumes that some extant species were not sampled.
}
\usage{
likelihood_bd(phylo, tot_time, f.lamb, f.mu, f, cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE,
expo.mu = FALSE, cond = "crown")
}

\arguments{
  \item{phylo}{  
phylo is an object of type phylogeny
}
  \item{tot_time}{
tot_time is the age of the phylogeny
}

  \item{f.lamb}{
The \eqn{\lambda} function. This corresponds to the speciation rate function at any time
The first argument is always time. The other arguments are the parameters of the variation,
given as a single vector.
}
  \item{f.mu}{
The \eqn{\mu} function. This corresponds to the extinction rate function at any time
The first argument is always time. The other arguments are the
parameters of the variation, given as a single vector.
}
  \item{f}{
Probability that some extant species were sampled and not included
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
Let \eqn{\Phi(s,t)} the probability that a lineage alive at time t has no descendant. Likelihood becomes:
\eqn{L(t_{1},...,t_{n})=\frac{f^n\Psi(t_{2}, t_{1})\prod_{i=2}^{n} \Psi(s_{i,1}, t_{i}) \Psi(s_{i,2}, t_{i}) }{1-\Phi(t_{1})}}{L(t[1],...,t[n])=f^n(\Psi(t[2], t[1])/(1-\Phi(t[1])))*\prod_{i=2}^{n} \Psi(s[i,1], t[i]) \Psi(s[i,2], t[i])}
where:
\itemize{
\item f denotes the probability that some extant species were sampled and not included in the sample
\item \eqn{s_{i,1},s_{i,2}}{s[i,1],s[i,2]} denote the times at which the daughter lineages introduced at time \eqn{t_{i}}{t[i]} themselves branch.}
}

\value{
  \item{loglikelihood}{The loglikelihood value}
}
\references{
Hélène Morlon, Todd L. Parsons, Joshua B. Plotkin:  \emph{Reconciling molecular phylogenies with the fossil record}, PNAS 2011

Fabien L. Condamine, Jonathan Rolland, Hélène Morlon: \emph{Macroevolutionary perspectives to environmental change}, Ecology Letters (2013)
}

\examples{
data(Cetacea)
tot_time <- max(node.age(Cetacea)$ages)
lamb_par <- c(1, -3)
# lambda(t) and mu(t)
f.lamb <- function(t){lamb_par[1] * exp(lamb_par[2] * t)}
f.mu <- function(t){0}
f <- 87/88
# model with none extinction & exponential speciation variation
likelihood <- likelihood_bd(Cetacea,tot_time,f.lamb,f.mu,f,cst.mu=TRUE,expo.lamb=TRUE)
print(likelihood)
# [1] -2160.751
}
