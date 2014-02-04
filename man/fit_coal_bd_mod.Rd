\name{fit_coal_bd_mod}
\alias{fit_coal_bd_mod}
\title{
Fit birth-death model using a coalescent approch
}
\description{
This function helps to fit various flavors of the birth-death model to a given phylogeny, by maximizing the likelihood function.
The default parameters allow to fit the most general model where the rates of speciation and extinction vary over time without a fixed extinction fraction and write as exponential models.
Outputs are the log-likelihood, the second order Akaike's Information Criterion, and the maximum likelihood estimates of the parameters of diversification. 
Depending on the model, these parameters include a combination of the speciation rate at present (lamb0), the exponential variation in speciation rate (alpha), the extinction rate at present (mu0), the exponential variation in extinction rate (beta) and the extinction fraction (extinction rate/speciation rate, eps). See notations in the PLoS Biology 2010 paper.
}

\usage{
fit_coal_bd_mod(phylo, lamb0 = 0.1, alpha = 1, mu0 = 0.01, beta = 0,
                meth = "Nelder-Mead", N0 = 0, cst.lamb = FALSE, cst.mu = FALSE,
                fix.eps = FALSE, mu.0 = FALSE, pos = TRUE, Vtimes = FALSE)
}

\arguments{
  \item{phylo}{
phylo is an object of type 'phylo' (see ape documentation)
}
  \item{lamb0}{
The speciation rate at present time.
}
  \item{alpha}{
The exponential variation in speciation rate. Model writes \eqn{\lambda(t) = lamb0 e^{alpha*t}}{lambda(t) = lamb0 exp(alpha t)}
}
  \item{mu0}{
The extinction rate at present time.
}
  \item{beta}{
The exponential variation in extinction rate. Model writes \eqn{mu(t) = mu0 e^{beta*t}}{mu(t) = mu0 exp(beta t)}
}
  \item{meth}{
The optimization method that should be used for the maximization of the likelihood function.
See \emph{optim} for more details.
}
  \item{N0}{
Number of species at present time
}
  \item{cst.lamb}{
Forces the speciation rate to be constant over time
}
  \item{cst.mu}{
Forces the extinction rate to be constant over time
}
  \item{fix.eps}{
Forces the extinction fraction to be constant over time
}
  \item{mu.0}{
Forces the extinction rate to 0
}
  \item{pos}{
Forces the rates of speciation and extinction to be positive.
Default value is TRUE
}
  \item{Vtimes}{
Time sorted from present to past
}

}
\details{
The function fits various flavors of the birth-death model to a given phylogeny, by maximizing the associated likelihood function.
The underlying birth-death models are described in PLoSB 2010 and may be resumed hereafter:
\itemize{\item Under model 3, \eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0}
\item Under model 4a, \eqn{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0}
\item Under model 4b, \eqn{\lambda(t) = \lambda_0, \mu(t) = \mu_0 exp(\beta t)}
\item Under model 4c, \eqn{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \epsilon \lambda(t)}
\item Under model 4d, \eqn{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = \mu_0 exp(\beta t)}
\item Under model 5, \eqn{\lambda(t) = \lambda_0, \mu(t) = 0}
\item Under model 6, \eqn{\lambda(t) = \lambda_0 exp(\alpha t), \mu(t) = 0}
}
where \eqn{\lambda_0, \alpha, \mu_0, \beta, \epsilon} are respectively:
\itemize{
\item The speciation rate at present;
\item The exponential variation in speciation rate;
\item The extinction rate at present;
\item The exponential variation in extinction rate;
\item The extinction fraction (extinction rate/speciation rate);
}
For more details, see notations in the PLoSB 2010 paper.

cst.lamb=TRUE forces the speciation rate to be constant over time (e.g. Models 3, 5 and 4b). 
cst.mu=TRUE forces the extinction rate to be constant over time (e.g. Models 3, and 4a). 
fix.eps forces the extinction fraction to be constant over time (e.g. Model 4c). 
mu.0=TRUE forces the extinction rate to 0 (e.g. Models 5 and 6).
pos=TRUE forces the rates of speciation and extinction to be positive.
}
\value{
  "fit_coal_bd_mod" returns a list containing the following components:
  \item{model}{The name of the used process}
  \item{LH}{The log-likelihood value}
  \item{aicc}{The second order Akaike's Information Criterion}
  \item{model.parameters}{Parameters of the models.}
}
\references{
Hélène Morlon, Matthew D. Potts, Joshua B. Plotkin: \emph{Inferring the Dynamics of Diversification: a Coalescent Approach}, PLoS Biology September 2010, Volume 8, Issue 9
}
\examples{
data(Cetaceae)
result <- fit_coal_bd_mod(Cetaceae, lamb0=0.01, alpha=-0.001, mu0=0.0, beta=0, N0=1e6)
print(result)
}
