\name{getMAPS_ClaDS}
\alias{getMAPS_ClaDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gets the Maximum A Posteriori for each ClaDS parameter
}
\description{
Extract the MAPs (Maximum A Posteriori) for the marginal posterior distributions estimated with fit_ClaDS
}
\usage{
getMAPS_ClaDS(sampler, burn = 1/2, thin = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampler}{
The output of a fit_ClaDS run.
}
  \item{burn}{
Number of iterations to drop in the beginning of the chains.}
  \item{thin}{
Thinning parameter, one iteration out of "thin" is kept to compute the MAPs.
}
}

\value{
A vector MAPS containing the MAPs for the marginal posterior distribution for each of the model's parameters.

MAPS[1:4] are the estimated hyperparameters, with MAPS[1] the sigma parameter (new rates stochasticity), MAPS[2] the alpha parameter (new rates trend),  MAPS[3] the turnover rate epsilon, and MAPS[4] the initial speciation rate lambda_0. 

MAPS[-(1:4)] are the estimated branch-specific speciation rates, given in the same order as the edges of the phylogeny on which the inference was performed.
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution}, doi 10.1038/s41559-019-0908-0
}
\author{
O. Maliet
}
\seealso{
  \code{\link{fit_ClaDS}}, \code{\link{plot_ClaDS_chains}}, \code{\link{getMAPS_ClaDS0}}
}
\examples{
data("Caprimulgidae_ClaDS2")

\donttest{
MAPS = getMAPS_ClaDS(Caprimulgidae_ClaDS2$sampler, thin = 1)

print(paste0("sigma = ", MAPS[1], " ; alpha = ", 
  MAPS[2], " ; epsilon = ", MAPS[3], " ; l_0 = ", MAPS[4] ))
plot_ClaDS_phylo(Caprimulgidae_ClaDS2$tree, MAPS[-(1:4)])
}
}
