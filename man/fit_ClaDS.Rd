\name{fit_ClaDS}
\alias{fit_ClaDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit ClaDS to a phylogeny
}
\description{
Performs the inferrence of branch-specific speciation rates and the model's hyper parameters for the model with constant extinction rate (ClaDS1) or constant turnover rate (ClaDS2).
}
\usage{
fit_ClaDS(tree,sample_fraction,iterations, thin = 50, file_name = NULL, it_save = 1000,
                     model_id="ClaDS2", nCPU = 1, mcmcSampler = NULL, ...)
}
\arguments{
  \item{tree}{
An object of class 'phylo'
}
  \item{sample_fraction}{
The sampling fraction for the clade on which the inference is performed.
}
  \item{iterations}{
Number of steps in the MCMC, should be a multiple of \code{thin}.
}
  \item{thin}{
Number of iterations between two chain state's recordings.
}
  \item{file_name}{
Name of the file in which the result will be saved. Use file_name = NULL (the default) to disable this option.}
  \item{it_save}{
Number of iteration between each backup of the result in file_name.
}
  \item{model_id}{
"ClaDS1" for constant extinction rate, "ClaDS2" (the default) for constant turnover rate.
}
  \item{nCPU}{
The number of CPUs to use. Should be either 1 or 3.
}
  \item{mcmcSampler}{
Optional output of \code{\link{fit_ClaDS}} to continue an already started run.}
  \item{\dots}{
Optional arguments, see details.}
}
\details{
This function uses a blocked Differential Evolution (DE) MCMC sampler, with sampling from the
past of the chains (Ter Braak, 2006; ter Braak and Vrugt, 2008). This sampler is self-adaptive
because proposals are generated from the past of the chains. In this sampler, three chains are
run simultaneously. Block updates is implemented by first drawing the number
of parameters to be updated from a truncated geometric distribution with mean 3, drawing
uniformly which parameter to update, and then following the normal DE algorithm.

The available optional arguments are : 
\describe{
  \item{Nchain}{
Number of MCMC chains (default to 3).
}
  \item{res_ClaDS0}{
The output of ClaDS0 to use as a startpoint. If NULL (the default) a random startpoint is used for the branch-specific speciation rates for each chain.
}
  \item{l0}{
The starting value for lambda_0 (not used if res_ClaDS0 != NULL).
}
  \item{s0}{
The starting value for sigma (not used if res_ClaDS0 != NULL).
}
  \item{nlambda}{
Number of subdivisions for the rate space discretization (use in the likelihood computation). Default to 1000.
}
  \item{nt}{
Number of subdivisions for the time space discretization (use in the likelihood computation). Default to 30.
}}
}

\value{
A 'list' object with fields :

\item{post }{The posterior function.}
\item{startvalue }{The starting value for the MCMC.}
\item{numPars }{The number of parameter in the model, including the branch-specific speciation rates.}
\item{Nchain }{The number of MCMC chains ran simultaneously.}
\item{currentLPs }{The current values of the logposterior for th \code{Nchains} chains.}
\item{proposalGenerator }{The proposal distribution for the MCMC sampler.}
\item{former }{The last output of \code{post} for each of the chains.}
\item{thin }{Number of iterations between two chain state's recordings.}
\item{alpha_effect }{A vector of size \code{nrow(tree$edge)}, where the ith element is the number of branches on the path from the crown of the tree and branch i (used internally in other functions).}
\item{consoleupdates }{The frequency at which the sampler state should be printed.}
\item{likelihood }{The likelihood function, used internally.}
\item{relToAbs }{A function mapping the relative changes in speciation rates to the absolute speciation rates for the object \code{phylo}, used internally.}
}
\references{
Ter Braak, C. J. 2006. A markov chain monte carlo version of the genetic algorithm differential
evolution: easy bayesian computing for real parameter spaces. \emph{Statistics and Computing 16:239-
249}.

ter Braak, C. J. and J. A. Vrugt. 2008. Differential evolution markov chain with snooker updater
and fewer chains. \emph{Statistics and Computing 18:435-446}.

Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution} in press
}
\author{
O. Maliet
}
\seealso{
\code{\link{fit_ClaDS0}}, \code{\link{plot_ClaDS_chains}}.
}
\examples{
\dontrun{
data("Caprimulgidae_ClaDS2")

tree = Caprimulgidae_ClaDS2$tree
sample_fraction = Caprimulgidae_ClaDS2$sample_fraction

sampler = fit_ClaDS(tree, sample_fraction, 1000, thin = 50, 
          file_name = NULL, model_id="ClaDS2", nCPU = 1)
plot_ClaDS_chains(sampler)

# continue the same run 
sampler = fit_ClaDS(tree, sample_fraction, 50, mcmcSampler = sampler)




# plot the result of the analysis (saved in "Caprimulgidae_ClaDS2")

data("Caprimulgidae_ClaDS2")

# plot the mcmc chains
plot_ClaDS_chains(Caprimulgidae_ClaDS2$sampler)

# extract the Maxima A Posteriori for each parameter
maps = getMAPS_ClaDS(Caprimulgidae_ClaDS2$sampler, thin = 1)
print(paste0("sigma = ", maps[1], " ; alpha = ", 
  maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] ))
  
# plot the infered branch specific speciation rates
plot_ClaDS_phylo(Caprimulgidae_ClaDS2$tree, maps[-(1:4)])
}

}
