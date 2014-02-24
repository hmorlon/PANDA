\name{likelihood_coal_var}
\alias{likelihood_coal_var}
\title{
Likelihood of a birth-death model using a coalescent approch
}
\description{
Computes the likelihood of a phylogeny under the 
Notations follow Morlon et al. PNAS 2011.
}
\usage{
likelihood_coal_var(Vtimes, ntips, lamb0, alpha, mu0, beta, N0, pos = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Vtimes}{
vector of times sorted from present to past
}
  \item{ntips}{
number of species in the phylogeny
}
  \item{lamb0}{
the speciation rate at present time
}
  \item{alpha}{
the exponential variation in speciation rate.  Model writes:
\deqn{\lambda(t) = \lambda_0 e^{\alpha t}}{\lambda(t) = \lambda_0 exp(\alpha t)}
}
  \item{mu0}{
extinction rate at present time
}
  \item{beta}{
exponential variation in speciation rate.  Model writes:
\deqn{\mu(t) = \mu_0 e^{\beta t}}{\mu(t) = \mu_0 exp(\beta t)}
}
  \item{N0}{
deterministic number of species at present time
}
  \item{pos}{
logical: should be set to FALSE to not enforce positive speciation and extinction rates
}
}

\value{
  a list containing the following components:
  \item{res}{the log-likelihood evaluation}
  \item{all}{the likelihood function evaluated for each time-length \eqn{t_i}{t[i]} between nodes i and i+1 of the clade}
}

\details{
Given a phylogeny, the function computes the likelihood of the distances \eqn{t_1,...,t_{k-2}}{t[1],...,t[k-2]} between
nodes in that phylogeny.
If we denote \eqn{t_i}{t[i]} the time-length between nodes i and i+1, and \eqn{v_i}{v[i]} the time-length between node i+1 and the present,
the likelihood writes as following:
\deqn{L(t_i) = i(i+1)\frac{\lambda(v_i)}{N(v_i)} exp(-i(i+1) \int_{v_i-t_i}^{v_i} \frac{\lambda(t)}{N(t)})}{L(t[i]) = i(i+1)(lambda(v[i]) / N(v[i])) exp(-i(i+1) integrate(\lambda(t)/N(t) dt, v[i]-t[i], v[i]) )}
The expression is valid if \eqn{\lambda(t) >= \mu(t)}.
With fixed speciation and extinction rate functions, log-likelihood is evaluated for each time-length \eqn{t_i}{t[i]} between nodes i and i+1 of the clade, and summed to give the total log-likelihood.

The likelihood function accounts \eqn{\lambda, \mu}, which might be of different types and may be resumed hereafter:
\itemize{
\item model 3: \eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0}{\lambda(t) = \lambda_0, \mu(t) = \mu_0}
(with cst.lamb=TRUE, cst.mu=TRUE)
\item model 4a: \eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \mu_0}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0}
(with cst.lamb=FALSE, cst.mu=TRUE)
\item model 4b: \eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0 e^{\beta t}}{\lambda(t) = \lambda_0, \mu(t) = \mu_0 exp(\beta t)}
(with cst.lamb=TRUE, cst.mu=FALSE)
\item model 4c: \eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \epsilon \lambda(t)}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \epsilon \lambda(t)}
(with cst.lamb=FALSE, cst.mu=FALSE, fix.eps=TRUE)
\item model 4d: \eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = \mu_0 e^{\beta t}}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0 exp(\beta t)}
(with cst.lamb=FALSE, cst.mu=FALSE, fix.eps=FALSE)
\item model 5: \eqn{\lambda(t) = \lambda_0, \mu(t) = 0}{\lambda(t) = \lambda_0, \mu(t) = 0}
(with cst.lamb=TRUE, mu.0=TRUE)
\item model 6: \eqn{\lambda(t) = \lambda_0 e^{\alpha t}, \mu(t) = 0}{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = 0}
(with cst.lamb=FALSE, mu.0=TRUE)
}
where \eqn{\lambda_0, \alpha, \mu_0, \beta, \epsilon} are respectively:
\itemize{
\item speciation rate at present;
\item exponential variation in speciation rate;
\item extinction rate at present;
\item exponential variation in extinction rate;
\item extinction fraction (extinction rate/speciation rate).
}

For more details, see notations in the PloSB 2010 paper.
}

\references{
Morlon, H., Potts, M.D., Plotkin, J.B. (2010) Inferring the Dynamics of Diversification: a Coalescent Approach, PLoS B 8(9): e1000493

Griffiths, R.C., Tavare, S.: Sampling theory for neutral alleles in a varying environment
}

\examples{
data(Cetacea)
Vtimes <- sort(branching.times(Cetacea))
# likelihood considering model 3
lamb0 <- 0.1
alpha <- -0.001
ntips <- Ntip(Cetacea)
likelihood <- likelihood_coal_var(Vtimes, ntips, lamb0, alpha, mu0=0, beta=0, N0=1e6)
}
