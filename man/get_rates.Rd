\name{get_rates}
\alias{get_rates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the absolute branch-specific rates.
}
\description{
Compute the branch-specific diversification rates from the relative changes in rates for a given phylogeny.}
\usage{
get_rates(phylo, lambda, log_scale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{
An object of class 'phylo'.
}
  \item{lambda}{
A vector of size \code{nrow(phylo$edge)+1}. The first element is lambda_0, the (i+1)th element is the relative change in rate for branch i (lambda_i / lambda_(parent_i)).}
  \item{log_scale}{
A boolean indicating whether \code{lambda} is given on a log scale or not. 
}
}
\value{
A vector of size \code{nrow(phylo$edge)+1}. The first element is lambda_0, the (i+1)th branch specific rate of branch i (lambda_i).
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution} in press
}
\author{
O. Maliet
}
\seealso{
\code{\link{sim_ClaDS}}, \code{\link{plot_ClaDS}}
}
\examples{
data(Cetacea)

relative_rates = rep(1, nrow(Cetacea$edge))
relative_rates[9] = 2
abs_rates = get_rates(Cetacea, c(0.1, relative_rates), log_scale = FALSE)
plot_ClaDS(Cetacea, abs_rates[-1])

abs_rates = get_rates(Cetacea, log(c(0.1, relative_rates)))
plot_ClaDS(Cetacea, exp(abs_rates[-1]))
}
