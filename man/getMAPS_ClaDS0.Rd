\name{getMAPS_ClaDS0}
\alias{getMAPS_ClaDS0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gets the Maximum A Posteriori for each ClaDS0 parameter
}
\description{
Extract the MAPs (Maximum A Posteriori) for the marginal posterior distributions estimated with run_ClaDS0.
}
\usage{
getMAPS_ClaDS0(phylo, sampler, burn=1/2, thin=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{
An object of class 'phylo'.
}
  \item{sampler}{
The output of a run_ClaDS0 run.
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
  \code{\link{fit_ClaDS0}}, \code{\link{plot_ClaDS0_chains}}, \code{\link{getMAPS_ClaDS}}
}
\examples{
set.seed(1)

obj= sim_ClaDS( lambda_0=0.1,    
                mu_0=0.5,      
                sigma_lamb=0.7,         
                alpha_lamb=0.90,     
                condition="taxa",    
                taxa_stop = 20,    
                prune_extinct = TRUE)  

tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]
data("ClaDS0_example")

# extract the Maximum A Posteriori for each of the parameters
MAPS = getMAPS_ClaDS0(ClaDS0_example$tree, 
                      ClaDS0_example$Cl0_chains, 
                      thin = 10)

# plot the simulated (on the left) and inferred speciation rates (on the right)
# on the same color scale
plot_ClaDS_phylo(ClaDS0_example$tree, 
          ClaDS0_example$speciation_rates, 
          MAPS[-(1:3)])
}
