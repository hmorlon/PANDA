\name{likelihood_coal_cst}
\alias{likelihood_coal_cst}

\title{
Likelihood of a constant diversity model with time-varying rate
}
\description{
Computes the likelihood of a phylogeny under assumption of constant diversity/exponentially turnover rate model.
Notations follow Morlon et al. PNAS 2011.
}
\usage{
likelihood_coal_cst(Vtimes, ntips, tau0, gamma, N0)
}
\arguments{
  \item{Vtimes}{
vector of times
}
  \item{ntips}{
number of tips
}
  \item{tau0}{
turnover rate at present time
}
  \item{gamma}{
exponential variation in turnover rate.
If gamma=0, the model is constant over time.
}
  \item{N0}{
number of species at present time
}
}

\value{
  a list containing the following components:
  \item{res}{the log-likelihood evaluation}
  \item{all}{the likelihood function evaluated for each time-length \eqn{t_i} between nodes i and i+1 of the clade}
}

\details{
Given a phylogeny, the function computes the likelihood of the distances \eqn{t_1,...,t_{k-2}} between
nodes in that phylogeny using a constant diversity hypothesis.
If we denote \eqn{t_i} the time-length between nodes i and i+1, \eqn{v_i} is the time-length between node i+1 and the present,
the likelihood writes as following:
%
\eqn{L(t_i) = i(i+1) \frac{\tau_0}{N_0} exp(-i(i+1) \frac{\tau_0}{N_0 \gamma} exp(\gamma v_i) (1-exp(-\gamma t_i)) ), \gamma \neq 0}
%
\eqn{L(t_i) = i(i+1) \frac{\tau_0}{N_0} exp(-i(i+1) \frac{\tau_0}{N_0} t_i), \gamma = 0}
%
where:
\itemize{
\item \eqn{\tau_0} is the turnover rate at present time;
\item \eqn{\gamma} is exponential variation in turnover rate.
}
}

\references{
Morlon, H., Potts, M.D., Plotkin, J.B. (2010) Inferring the Dynamics of Diversification: a Coalescent Approach, PLoS B 8(9): e1000493

Griffiths, R.C., Tavare, S.: Sampling theory for neutral alleles in a varying environment
}

\examples{
data(Cetacea)
Vtimes <- sort(branching.times(Cetacea))
tau0 <- 0.1
gamma <- -0.001
ntips <- Ntip(Cetacea)
N0 <- 1e6
likelihood <- likelihood_coal_cst(Vtimes,ntips,tau0, gamma,N0)
}
