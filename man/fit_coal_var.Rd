\name{fit_coal_var}
\alias{fit_coal_var}
\title{
Fit birth-death model using a coalescent approch
}
\description{
Fits various flavors of the birth-death model to a phylogeny, by maximum likelihood. Notations follow Morlon et al. PNAS 2011.
}

\usage{
fit_coal_var(phylo, lamb0 = 0.1, alpha = 1, mu0 = 0.01, beta = 0,
             meth = "Nelder-Mead", N0 = 0, cst.lamb = FALSE, cst.mu = FALSE,
             fix.eps = FALSE, mu.0 = FALSE, pos = TRUE, Vtimes = FALSE)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{lamb0}{
the speciation rate at present time
}
  \item{alpha}{
the exponential variation in speciation rate. Model writes \eqn{\lambda(t) = lamb0 e^{\alpha*t}}{lambda(t) = lamb0 exp(alpha t)}
}
  \item{mu0}{
the extinction rate at present time
}
  \item{beta}{
the exponential variation in extinction rate. Model writes \eqn{\mu(t) = mu0 e^{\beta*t}}{mu(t) = mu0 exp(beta t)}
}
  \item{meth}{
optimization to use to maximize the likelihood function, see \emph{optim} for more details.
}
  \item{N0}{
number of species at present time
}
  \item{cst.lamb}{
logical: should be set to TRUE only if f.lamb is constant (i.e. does not depend on time) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{cst.mu}{
logical: should be set to TRUE only if f.mu is constant (i.e. does not depend on time) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{fix.eps}{
logical: should be set to TRUE only if extinction fraction is constant (i.e. does not depend on time)
}
  \item{mu.0}{
logical: should be set to TRUE to force the extinction rate to 0
}
  \item{pos}{
logical: should be set to FALSE to not enforce positive speciation and extinction rates
}
  \item{Vtimes}{
time sorted from present to past
}
}

\value{
	a list with the following components
  \item{model}{The name of the used process}
  \item{LH}{The log-likelihood value}
  \item{aicc}{The second order Akaike's Information Criterion}
  \item{model.parameters}{Parameters of the model}
}

\details{
The function fits various flavors of the birth-death model to a given phylogeny, by maximizing the associated likelihood function.
The default parameters allow to fit the most general model where speciation and extinction rates vary over time without a fixed extinction fraction and write as exponential models.
Depending on the model, these parameters include a combination of the speciation rate at present (lamb0), the exponential variation in speciation rate (alpha), the extinction rate at present (mu0), the exponential variation in extinction rate (beta) and the extinction fraction (extinction rate/speciation rate, eps). See notations in the PLoS Biology 2010 paper.

The underlying birth-death models may be resumed hereafter:
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
}

\references{
Morlon, H., Potts, M.D., Plotkin, J.B. (2010) Inferring the Dynamics of Diversification: a Coalescent Approach, PLoS B 8(9): e1000493
}

\seealso{
\code{\link[PANDA]{likelihood_coal_var}}
}

\examples{
data(Cetacea)
result <- fit_coal_var(Cetacea, lamb0=0.01, alpha=-0.001, mu0=0.0, beta=0, N0=1e6)
print(result)
}
