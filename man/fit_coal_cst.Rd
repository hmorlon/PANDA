\name{fit_coal_cst}
\alias{fit_coal_cst}

\title{
Fit diversity model with time-varying rate
}
\description{
Fits the constant diversity/exponentially turnover rate model to a phylogeny, by maximum likelihood. Notations follow Morlon et al. PLoSB 2010.
}
\usage{
fit_coal_cst(phylo, tau0 = 1e-2, gamma = 1, cst.rate = FALSE,
             meth = "Nelder-Mead", N0 = 0, Vtimes = FALSE)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{tau0}{
initial turnover rate at present time
}
  \item{gamma}{
initial exponentially variation in turnover rate
}
  \item{cst.rate}{
logical: should be set to TRUE to fit model with time-constant rate.  By default, a model with time-varying rate (exponential) is used.
}
  \item{meth}{
optimization to use to maximize the likelihood function, see \emph{optim} for more details.
}
  \item{N0}{
initial number of species.
With default value(0), the parameter is set to the number of tips in the phylogeny (i.e. the phylogeny is assumed to be 100\% complete)
}
  \item{Vtimes}{
logical: should be set to TRUE to take branching times into account.  Otherwise, the phylogeny is taken into account.
}
}

\value{
	a list with the following components
  \item{model}{The name of the used process}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{tau0}{the estimated turnover rate at present time}
  \item{gamma}{exponential variation in turnover rate (if cst.rate is FALSE)}
}

\details{
The function fits the constant diversity/exponentially turnover rate model to a given phylogeny, by maximizing the associated likelihood function.
Assuming we know the total number of species at present, we estimate the turnover rate parameters. 
The number of species is set to the number of tips in the phylogeny if default value (0) is fixed, assuming
the phylogeny be 100\% complete.
In addition to the turnover rate at present/the exponential variation in turnover rate, the log-likelihood and the second order Akaike's Information Criterion are returned.
The underlying model (Hey/Moran Model) may be resumed hereafter:
\eqn{ \tau(t) = \tau_0 e^{-\gamma t}}{ \tau(t) = \tau_0 exp(-\gamma t)}
where \eqn{\tau(t)} is the turnover rate function.
Under these assumptions, the model is said to be satured diversity with positive extinction.
For more details, see notations in the PloSB 2010 paper (model 2).

The meth argument fixes the optimization method. Depending on this choice, one should take
care of the default values of tau0 (1e-2) and gamma (1). Indeed, the likelihood function has several local maxima.
To make sure to not be stuck in a local optimum, several inital values might be investigated.
}

\references{
Morlon, H., Potts, M.D., Plotkin, J.B. (2010) Inferring the Dynamics of Diversification: a Coalescent Approach, PLoS B 8(9): e1000493
}

\seealso{
\code{\link[PANDA]{likelihood_coal_cst}}
}

\examples{
data(Cetacea)
result <- fit_coal_cst(Cetacea, tau0=1.e-3, gamma=-1, cst.rate=FALSE)
print(result)
}
