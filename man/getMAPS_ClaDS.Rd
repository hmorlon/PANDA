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
A vector containing the MAPs for the marginal posterior distribution for each of the model's parameters.
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution} in press
}
\author{
O. Maliet
}
\seealso{
  \code{\link{fit_ClaDS}}, \code{\link{plot_ClaDS_chains}}, \code{\link{getMAPS_ClaDS0}}
}
\examples{
data("Caprimulgidae_ClaDS2")

maps = getMAPS_ClaDS(Caprimulgidae_ClaDS2$sampler, thin = 5)

print(paste0("sigma = ", maps[1], " ; alpha = ", 
  maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] ))
plot_ClaDS_phylo(Caprimulgidae_ClaDS2$tree, maps[-(1:4)])

}
