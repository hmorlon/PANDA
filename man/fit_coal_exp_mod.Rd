\name{fit_coal_exp_mod}
\alias{fit_coal_exp_mod}

\title{
Fit diversity model with time-varying rate
}
\description{
This function helps to fit the constant diversity/exponentially turnover rate model to a given phylogeny, by maximizing the likelihood function.
Assuming we know the total number of species at present, we estimate the turnover rate parameters. 
The number of species is set to the number of tips in the phylogeny if default value (0) is fixed, assuming
the phylogeny be 100\% complete.
In addition to the turnover rate at present/the exponential variation in turnover rate, the log-likelihood and the second order Akaike's Information Criterionare returned.
}
\usage{
fit_coal_exp_mod(phylo, tau0 = 1e-2, gamma = 1, meth = "Nelder-Mead", N0 = 0, Vtimes = FALSE)
}

\arguments{
  \item{phylo}{
phylo is an object of type phylogeny
}
  \item{tau0}{
Initial turnover rate at present time
}
  \item{gamma}{
Initial exponentially variation in turnover rate
}
  \item{meth}{
The optimization method that should be used for the maximization of the likelihood function.
See \emph{optimize} for more details.
}
  \item{N0}{
Initial number of species.
With default value(0), the parameter is set to the number of tips in the phylogeny (i.e. the phylogeny is assumed to be 100% complete)
}
  \item{Vtimes}{
  Boolean. If TRUE, the function takes into account branching times. Otherwise, the phylogeny is taken into account.
}
}
\details{
The function fits the constant diversity/exponentially turnover rate model to a given phylogeny, by maximizing the associated likelihood function.
The underlying model (Hey/Moran Model) may be resumed hereafter:
\eqn{ \tau(t) = \tau_0 exp(-\gamma t)}
where \eqn{\tau(t)} is the turnover rate function.
Under these assumtions, the model is said to be satured diversity with positive extinction.
For more details, see notations in the PloSB 2010 paper (model 2).
The meth argument fixes the optimization method. Depending on the choose of this one, one should
care with the default initial value of tau0 (1e-2) and gamma (1). Indeed, the likelihood function has several local maxima.
To make sure to not be stuck in a local optimum, several inital values might be investigated.
The number of species is set to the number of tips in the phylogeny if default value (0) is fixed, assuming
the phylogeny is 100\% complete.
}
\value{
  "fit_coal_exp_mod" returns a list containing the following components:
  \item{model}{The name of the used process}
  \item{LH}{The log-likelihood value}
  \item{aicc}{The second order Akaike's Information Criterion}
  \item{tau0}{The estimated turnover rate at present time}
  \item{gamma}{Exponential variation in turnover rate}
}
\references{
H.Morlon, M.D.Plotts, J.B.Plotkin: Inferring the Dynamics of Diversification: a Coalescent Approach, PloSB 2010
}
\author{
H.Morlon, \email{helene.morlon@cmap.polytechnique.fr}
}
\examples{
data(Cetaceae)
result <- fit_coal_exp_mod(Cetaceae, tau0=1.e-3, gamma=-1)
print(result)
}
