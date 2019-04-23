\name{prepare_ClaDS}
\alias{prepare_ClaDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare the inference of ClaDS parameters}
\description{
Prepare the inference for the model with constant extinction rate (ClaDS1) or constant turnover rate (ClaDS2).
}
\usage{
prepare_ClaDS(tree, sample_fraction, Nchain = 3, 
                model_id = "ClaDS2", res_ClaDS0 = NULL, l0 = 0.1, 
                s0 = 1, nlambda = 1000, nt = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
An object of class 'phylo'.
}
  \item{sample_fraction}{
The sampling fraction for the clade on which the inference is performed.
}
  \item{Nchain}{
Number of MCMC chains (default to 3).
}
  \item{model_id}{
"ClaDS1" for constant extinction rate, "ClaDS2" (the default) for constant turnover rate.
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
}
}
\value{
A 'list' object with fields :
\item{likelihood }{The likelihood function as returned by the function \code{\link{createLikelihood_ClaDS2}}.}
\item{relToAbs }{A function mapping the relative changes in speciation rates to the absolute speciation rates for the object \code{phylo}, as returned by the function \code{\link{createLikelihood_ClaDS2}}.}
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
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution} in press
}
\author{
O. Maliet
}
\seealso{
  \code{\link{run_ClaDS}}, \code{\link{add_iter_ClaDS}}, \code{\link{plotChains_ClaDS}}, \code{\link{getMAPS_ClaDS}}
}
\examples{

data("Caprimulgidae_ClaDS2")

sampler = prepare_ClaDS(Caprimulgidae_ClaDS2$tree,Caprimulgidae_ClaDS2$sample_fraction,l0=0.1)
}
