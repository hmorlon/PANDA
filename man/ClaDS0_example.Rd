\name{ClaDS0_example}
\alias{ClaDS0_example}
\docType{data}
\title{
An example run of ClaDS0.
}
\description{
An example of the run on the inference of ClaDS0 on a simulated phylogeny, thinned every 10 iterations.}
\usage{data("ClaDS0_example")}
\format{
A list object with fields :

\describe{
\item{\code{tree}}{The simulated phylogeny on which we ran the model.}
\item{\code{speciation_rates}}{The simulated speciation rates.}
\item{\code{Cl0_chains}}{The output of the \code{run_ClaDS0} run.}
}
}
\seealso{\code{\link{fit_ClaDS0}}}

\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution} in press
}
\examples{
data(ClaDS0_example)

# plot the resulting chains for the first 4 parameters
plot_ClaDS0_chains(ClaDS0_example$Cl0_chains, param = 1:4)

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
