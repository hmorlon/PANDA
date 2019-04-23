\name{createLikelihood_ClaDS2}
\alias{createLikelihood_ClaDS2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood of a tree for the ClaDS2 model
}
\description{
Prepares the computation of the likelihood of a phylogeny for the model with constant turnover rate.
}
\usage{
createLikelihood_ClaDS2(phylo, root_depth = 0, relative = F, nt = 1000,
                    nlambda = 100, mlambda = 1e-25, Mlambda = 1e+10, conv
                    = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{
An object of class 'phylo'.
}
  \item{root_depth}{
Time before the first speciation event in the reconstructed phylogeny (default to 0, in which case the likelihood is conditionned with the crown age). }
  \item{relative}{
A boolean indicating whether the speciation rates are given as relative changes in rates (TRUE, log(lambda/lambda_parent)) or as the absolute rate value lambda (FALSE, the default). }
  \item{nt}{
Number of subdivisions for the time space discretization (use in the likelihood computation). Default to 30.
}
  \item{nlambda}{
Number of subdivisions for the rate space discretization (use in the likelihood computation). Default to 1000.
}
  \item{mlambda}{
Minimum possible value of the speciation rates.
}
  \item{Mlambda}{
Maximum possible value of the speciation rates.
}
  \item{conv}{
The required relative precision in the Magnus expansion.
}
}

\value{
A 'list' object with two elements
\item{ll }{A function with two arguments, \code{param} (the model parameters) and \code{f} that compute the loglikelihood of \code{phylo} for the ClaDS1 model}
\item{relToAbs }{A function mapping the relative changes in speciation rates to the absolute speciation rates for the object \code{phylo} (see the \code{\link{get_rates}} function).}
}
\note{
The likelihood is computed through the discretization of the rate space and the use of a Magnus expansion (Maliet et al. 2019). 
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution} in press
}
\author{
O. Maliet
}
\seealso{
   \code{\link{createLikelihood_ClaDS2}}, \code{\link{prepare_ClaDS}}, \code{\link{run_ClaDS}}, \code{\link{add_iter_ClaDS}}
}
\examples{
data("Cetacea")
param = c(0.5,1,0.8,rep(0.1,nrow(Cetacea$edge)+1))

fun = createLikelihood_ClaDS2(Cetacea)
print(fun$ll(param,1)$LL)
}