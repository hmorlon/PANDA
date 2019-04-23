\name{add_iter_ClaDS}
\alias{add_iter_ClaDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add iterations to a ClaDS fit
}
\description{
Add iterations to the inferrence of branch-specific speciation rates and the model's hyper parameters for the model with constant extinction rate (ClaDS1) or constant turnover rate (ClaDS2).
}
\usage{
add_iter_ClaDS(mcmcSampler, iterations, thin = NULL, nCPU = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcSampler}{
The output of prepare_ClaDS or fit_ClaDS
}
  \item{iterations}{
Number of steps in the MCMC.
}
  \item{thin}{
Number of iterations between two chain state's recordings.
}
  \item{nCPU}{
The number of CPUs to use. Should be either 1 or 3.
}
}
\details{
This function is used internally within the \code{run_ClaDS} function.

It uses a blocked Differential Evolution (DE) MCMC sampler, with sampling from the
past of the chains (Ter Braak, 2006; ter Braak and Vrugt, 2008). This sampler is self-adaptive
because proposals are generated from the past of the chains. In this sampler, three chains are
run simultaneously. Block updates is implemented by first drawing the number
of parameters to be updated from a truncated geometric distribution with mean 3, drawing
uniformly which parameter to update, and then following the normal DE algorithm.}

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
Ter Braak, C. J. 2006. A markov chain monte carlo version of the genetic algorithm dierential
evolution: easy bayesian computing for real parameter spaces. \emph{Statistics and Computing 16:239-
249}.

ter Braak, C. J. and J. A. Vrugt. 2008. Dierential evolution markov chain with snooker updater
and fewer chains. \emph{Statistics and Computing 18:435-446}.

Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution} in press
}
\author{
O. Maliet
}

\seealso{
\code{\link{prepare_ClaDS}}, \code{\link{run_ClaDS}}.
}
\examples{
\dontrun{
data("Caprimulgidae_ClaDS2")

tree = Caprimulgidae_ClaDS2$tree
sample_fraction = Caprimulgidae_ClaDS2$sample_fraction

sampler = add_iter_ClaDS(Caprimulgidae_ClaDS2$sampler, 100)
}
}
