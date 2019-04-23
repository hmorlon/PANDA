\name{plotChains_ClaDS0}
\alias{plotChains_ClaDS0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the MCMC chains obtained when infering ClaDS0 parameters
}
\description{
Plot the MCMC chains obtained with run_ClaDS0.
}
\usage{
plotChains_ClaDS0(sampler, burn = 1/2, thin = 1, param = c("sigma", "alpha", "l_0", "LP"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{sampler}{
The output of a run_ClaDS0 run.
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
  \code{\link{run_ClaDS0}}, \code{\link{getMAPS_ClaDS0}}, \code{\link{plotChains_ClaDS}}
}
\examples{
data("ClaDS0_example")

plotChains_ClaDS0(ClaDS0_example$Cl0_chains)
plotChains_ClaDS0(ClaDS0_example$Cl0_chains, param = paste0("lambda_", c(1,10,5)))

}
