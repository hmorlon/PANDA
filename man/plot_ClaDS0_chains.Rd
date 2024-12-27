\name{plot_ClaDS0_chains}
\alias{plot_ClaDS0_chains}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the MCMC chains obtained when infering ClaDS0 parameters
}
\description{
Plot the MCMC chains obtained with run_ClaDS0.
}
\usage{
plot_ClaDS0_chains(sampler, burn = 1/2, thin = 1, 
                  param = c("sigma", "alpha", "l_0", "LP"))
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

\value{
Plot representing parameter MCMC chains
}

\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution}, doi 10.1038/s41559-019-0908-0
}
\author{
O. Maliet
}
\seealso{
  \code{\link{fit_ClaDS0}}, \code{\link{getMAPS_ClaDS0}}, \code{\link{plot_ClaDS_chains}}
}
\examples{
data("ClaDS0_example")

plot_ClaDS0_chains(ClaDS0_example$Cl0_chains)
plot_ClaDS0_chains(ClaDS0_example$Cl0_chains, param = paste0("lambda_", c(1,10,5)))

}
