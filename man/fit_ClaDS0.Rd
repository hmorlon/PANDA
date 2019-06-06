\name{fit_ClaDS0}
\alias{fit_ClaDS0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer ClaDS0's parameter on a phylogeny}
\description{
Infer branch-specific speciation rates and the model's hyper parameters for the pure-birth model
}
\usage{
fit_ClaDS0(tree, name, pamhLocalName = "pamhLocal",  
            iteration = 1e+07, thin = 20000, update = 1000, 
            adaptation = 10, seed = NULL, nCPU = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
An object of class 'phylo'.
}
  \item{name}{
The name of the file in which the results will be saved. Use name = NULL to disable this option.
}
  \item{pamhLocalName}{
The function is writing in a text file to make the execution quicker, this is the name of this file.
}
  \item{iteration}{
Number of iterations after which the gelman factor is computed and printed. The function stops if it is below 1.05
}
  \item{thin}{
Number of iterations between two chain state's recordings.
}
  \item{update}{
Number of iterations between two adjustments of the proposal parameters during the adaptation phase of the sampler.
}
  \item{adaptation}{
Number of times the proposal is adjusted during the adaptation phase of the sampler.
}
  \item{seed}{
An optional seed for the MCMC run.
}
  \item{nCPU}{
The number of CPUs to use. Should be either 1 or 3.
}
}
\details{
This function uses a Metropolis within Gibbs MCMC sampler with a bactrian proposal (ref) with an initial adaptation phase. During this phase, the proposal is adjusted "adaptation" times every "update" iterations to reach a goal acceptance rate of 0.3.

To monitor convergence, 3 independant MCMC chains are run simultaneously and the Gelman statistics is computed every "iteration" iterations. The inference is stopped when the maximum of the one dimentional Gelman statistics (computed for each of the parameters) is below 1.05.}
\value{
A mcmc.list object with the three MCMC chains.
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution}, doi 10.1038/s41559-019-0908-0

}
\author{
O. Maliet
}
\seealso{
   \code{\link{getMAPS_ClaDS0}}, \code{\link{plot_ClaDS0_chains}}, \code{\link{fit_ClaDS}}
}
\examples{
\donttest{

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

plot_ClaDS_phylo(tree,speciation_rates)

sampler = fit_ClaDS0(tree=tree,        
              name="ClaDS0_example.Rdata",      
              nCPU=1,               
              pamhLocalName = "local",
              iteration=500000,
              thin=2000,
              update=1000, adaptation=5) 
              
}}
