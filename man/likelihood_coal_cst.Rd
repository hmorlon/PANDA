\name{likelihood_coal_cst}
\alias{likelihood_coal_cst}

\title{
Likelihood of a constant diversity model with time-varying rate
}
\description{
This function computes the log-likelihood to a given phylogeny under assumption of constant diversity/exponentially turnover rate model.
With a fixed number of species at present time, the log-likelihood is evaluated for each time-length \eqn{t_i} between nodes i and i+1 of the clade
then final log-likelihood is given as sum of these last ones.
}
\usage{
likelihood_coal_cst(Vtimes, ntips, tau0, gamma, N0)
}
\arguments{
  \item{Vtimes}{
Vector of times
}
  \item{ntips}{
Number of tips
}
  \item{tau0}{
Turnover rate at present time
}
  \item{gamma}{
Exponential variation in turnover rate
If gamma=0, the model is time-constant rate
}
  \item{N0}{
Number of species at present time
}
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
\item \eqn{\gamma} is exponential variation in turnover rate;
}
}
\value{
  "likelihood_coal_cst" returns a list containing the following components:
  \item{res}{The log-likelihood evaluation}
  \item{all}{The likelihood function evaluated for each time-length \eqn{t_i} between nodes i and i+1 of the clade}
}
\references{
H.Morlon, M.D.Plotts, J.B.Plotkin: Inferring the Dynamics of Diversification: a Coalescent Approach, PloSB 2010
R.C.Griffiths, S.Tavare: Sampling theory for neutral alleles in a varying environment
}
\author{
H.Morlon, \email{helene.morlon@cmap.polytechnique.fr}
}
\examples{
data(bird.families)
Vtimes <- sort(branching.times(bird.families))
tau0 <- 0.1
gamma <- -0.001
ntips <- Ntip(bird.families)
N0 <- 1e6
likelihood <- likelihood_coal_cst(Vtimes,ntips,tau0, gamma,N0)
}
