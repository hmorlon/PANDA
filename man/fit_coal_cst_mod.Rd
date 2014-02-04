\name{fit_coal_cst_mod}
\alias{fit_coal_cst_mod}

\title{
Fit diversity model with time-constant rate
}
\description{
This function helps to fit the constant diversity/constant turnover rate model to a given phylogeny, by maximizing the likelihood function.
Assuming we know the total number of species at present, we estimate the turnover rate. 
The number of species is set to the number of tips in the phylogeny if default value (0) is fixed, assuming
the phylogeny is 100\% complete.
In addition to the turnover rate parameter, the log-likelihood and the second order Akaike's Information Criterion are returned.
}
\usage{
fit_coal_cst_mod(phylo, tau0 = 1, meth = "Nelder-Mead", N0 = 0, Vtimes = FALSE)
}

\arguments{
  \item{phylo}{
phylo is an object of type phylogeny
}
  \item{tau0}{
Initial turnover rate. The fit method will estimate the one which maximizes
the likelihood function.
}
  \item{meth}{
The optimization method that should be used for the maximization of the likelihood function.
See \emph{optimize} for more details.
}
  \item{N0}{
Initial number of species.
With default value(0), the parameter is set to the number of tips in the phylogeny (i.e. the phylogeny is assumed to be 100\% complete)
}
  \item{Vtimes}{
Boolean. If TRUE, the function takes into account branching times. Otherwise, the phylogeny is taken into account.
}
}
\details{
The function fits the constant diversity/constant turnover rate model to a given phylogeny, by maximizing the associated likelihood function.
The underlying model (Hey/Moran Model) may be resumed hereafter:
\eqn{ \tau(t) = \tau_0}
where \eqn{\tau(t)} is the turnover rate function.
Under these assumtions, the model is said to be satured diversity with positive extinction.
For more details, see notations in the PLoS Biology 2010 paper.
The meth argument fixes the optimization method. Depending on the choice of this parameter, one should
care with the default initial value of tau0 (1). Indeed, the likelihood function has several local maxima.
To make sure to not be stuck in a local optimum, several inital values might be investigated.
The number of species is set to the number of tips in the phylogeny if default value (0) is fixed, assuming
the phylogeny is 100\% complete.
}
\value{
  "fit_coal_cst_mod" returns a list containing the following components:
  \item{model}{The name of the used process}
  \item{LH}{The log-likelihood value}
  \item{aicc}{The second order Akaike's Information Criterion}
  \item{tau0}{The estimated turnover rate at present time}
}
\references{
Hélène Morlon, Matthew D. Potts, Joshua B. Plotkin: \emph{Inferring the Dynamics of Diversification: a Coalescent Approach}, PLoS Biology September 2010, Volume 8, Issue 9
}

\examples{
data(Cetaceae)
result <- fit_coal_cst_mod(Cetaceae, tau0=1.5, N0=1e6)
print(result)
}
