\name{plotChains_ClaDS}
\alias{plotChains_ClaDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the MCMC chains obtained when infering ClaDS parameters}
\description{
Plot the MCMC chains obtained with fit_ClaDS.
}
\usage{
plotChains_ClaDS(sampler, burn = 1/2, thin = 1, param = c("sigma", "alpha", "mu", "LP"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampler}{
The output of a fit_ClaDS run.
}
  \item{burn}{
Number of iterations to drop in the beginning of the chains.
}
  \item{thin}{
Thinning parameter, one iteration out of "thin" is plotted.
}
  \item{param}{
Either a vector of "character" elements with the name of the parameter to plot, or a vector of integers indicating what parameters to plot. 
}
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution} in press
}
\author{
O. Maliet
}
\seealso{
  \code{\link{run_ClaDS}}, \code{\link{getMAPS_ClaDS}}, \code{\link{plotChains_ClaDS0}}
}
\examples{
data("Caprimulgidae_ClaDS2")

plotChains_ClaDS(Caprimulgidae_ClaDS2$sampler)

plotChains_ClaDS(Caprimulgidae_ClaDS2$sampler, burn = 1/4, param = c("sigma", "alpha", "l_0", "LP"))

plotChains_ClaDS(Caprimulgidae_ClaDS2$sampler, burn = 1/5, thin = 5, param = c(1,5,6,15))

}
