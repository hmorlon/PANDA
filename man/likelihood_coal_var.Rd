\name{likelihood_coal_var}
\alias{likelihood_coal_var}
\title{
Likelihood of a birth-death model using a coalescent approch
}
\description{
This function computes the log-likelihood of a given phylogeny 
under assumption of positive ratio (\eqn{\lambda(t) \geq \mu(t)}) and various kinds of rate variation \eqn{\lambda, \mu}.
With fixed speciation and extinction rate functions, log-likelihood is evaluated for each time-length \eqn{t_i}{t[i]} between nodes i and i+1 of the clade, and summed to give the total log-likelihood.
}
\usage{
likelihood_coal_var(Vtimes, ntips, lamb0, alpha, mu0, beta, N0, pos = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Vtimes}{
Vector of times sorted from present to past
}
  \item{ntips}{
Number of species in the considered phylogeny
}
  \item{lamb0}{
The speciation rate at present time
}
  \item{alpha}{
The exponential variation in speciation rate. 
Model writes:
\deqn{\lambda(t) = \lambda_0 exp(\alpha t)}
}
  \item{mu0}{
The extinction rate at present time
}
  \item{beta}{
The exponential variation in speciation rate. 
Model writes:
\deqn{\mu(t) = \mu_0 exp(\beta t)}
}
  \item{N0}{
Deterministic number of species at present time
}
  \item{pos}{
Boolean, forces the rates of speciation and extinction to be positive (TRUE)
}
}
\details{
Given a phylogeny, the function computes the likelihood of the distances \eqn{t_1,...,t_{k-2}}{t[1],...,t[k-2]} between
nodes in that phylogeny.
If we denote \eqn{t_i}{t[i]} the time-length between nodes i and i+1, and \eqn{v_i}{v[i]} the time-length between node i+1 and the present,
the likelihood writes as following:
\deqn{L(t_i) = i(i+1)\frac{\lambda(v_i)}{N(v_i)} exp(-i(i+1) \int_{v_i-t_i}^{v_i} \frac{\lambda(t)}{N(t)})}{L(t[i]) = i(i+1)(lambda(v[i]) / N(v[i])) exp(-i(i+1) integrate(\lambda(t)/N(t) dt, v[i]-t[i], v[i]) )}
The expression is valid if \eqn{\lambda(t) >= \mu(t)}.
The likelihood function accounts \eqn{\lambda, \mu}, which might be of different types and may be resumed hereafter:
\itemize{\item Under model 3, \eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0}
\item Under model 4a, \eqn{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0}
\item Under model 4b, \eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0 exp(\beta t)}
\item Under model 4c, \eqn{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \epsilon \lambda(t)}
\item Under model 4d, \eqn{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0 exp(\beta t)}
\item Under model 5, \eqn{\lambda(t) = \lambda_0, \mu(t) = 0}
\item Under model 6, \eqn{\lambda(t) = \lambda_0  exp(\alpha t), \mu(t) = 0}
}
where \eqn{\lambda_0, \alpha, \mu_0, \beta, \epsilon} are respectively:
\itemize{
\item The speciation rate at present;
\item The exponential variation in speciation rate;
\item The extinction rate at present;
\item The exponential variation in extinction rate;
\item The extinction fraction (extinction rate/speciation rate);
}
For more details, see notations in the PloSB 2010 paper.

cst.lamb=TRUE forces the speciation rate to be constant over time (e.g. Models 3, 5 and 4b). 
cst.mu=TRUE forces the extinction rate to be constant over time (e.g. Models 3, and 4a). 
fix.eps forces the extinction fraction to be constant over time (e.g. Model 4c). 
mu.0=TRUE forces the extinction rate to 0 (e.g. Models 5 and 6).
pos=TRUE forces the rates of speciation and extinction to be positive over the history of the clade.
}
\value{
  "likelihood_coal_var" returns a list containing the following components:
  \item{res}{The log-likelihood evaluation}
  \item{all}{The likelihood function evaluated for each time-length \eqn{t_i}{t[i]} between nodes i and i+1 of the clade}
}
\references{
Hélène Morlon, Matthew D. Potts, Joshua B. Plotkin: \emph{Inferring the Dynamics of Diversification: a Coalescent Approach}, PLoS Biology September 2010, Volume 8, Issue 9

R. C. Griffiths, S. Tavare: \emph{Sampling theory for neutral alleles in a varying environment}
}

\examples{
data(bird.families)
Vtimes <- sort(branching.times(bird.families))
# likelihood considering model 3
lamb0 <- 0.1
alpha <- -0.001
ntips <- Ntip(bird.families)
likelihood <- likelihood_coal_var(Vtimes, ntips, lamb0, alpha, mu0=0, beta=0, N0=1e6)
}
