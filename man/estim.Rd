\name{estim}
\alias{estim}

\title{
Estimation of traits ancestral states.
}
\description{
Reconstruct the ancestral states at the root and for each nodes of phylogenetic tree from parameters obtained using the \code{fit_t_pl} function.}
\usage{

estim(object)
  
}
\arguments{
  \item{object}{
A penalized likelihood model fit obtained by the \code{fit_t_pl} function.
}
 

\value{
	a list with the following components
  
  \item{root}{the reconstructed ancestral states at the root}
  \item{nodes}{the reconstructed ancestral states at each nodes}

  }

\details{
\code{estim} reconstruct the ancestral states at the root and for each nodes of a phylogenetic tree from penalized likelihood models fit to high-dimensional datasets (where the number of variables p is larger than n; see details for the models options in \code{fit_t_pl}). Ancestral states are estimated using generalized least squares (GLS; Martins & Hansen 1997, Cunningham et al. 1998 ).}

\note{
	Contrary to conventional PCA, the principal axes of the phylogenetic PCA are not orthogonal, they represent the main axes of (independent) evolutionary changes.}
\references{
Clavel, J. Aristide, L., Morlon, H., 2017. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. submitted

Cunningham C.W., Omland K.E., Oakley T.H. 1998. Reconstructing ancestral character states: a critical reappraisal. Trends Ecol. Evol. 13:361-366.

Martins E.P., Hansen T.F. 1997. Phylogenies and the comparative method: a general approach to incorporating phylogenetic information into the analysis of interspecific data. Am. Nat. 149:646-667.
}

\author{J. Clavel}

\seealso{
\code{\link{fit_t_pl}},
\code{\link{penalized_phyl.pca}},
\code{\link{GIC}}, 
\code{\link{gic_criterion}}
}

\examples{
\donttest{
require(mvMORPH)
set.seed(1)
n <- 32 # number of species
p <- 31 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- Posdef(p)      # a random symmetric matrix (covariance)

# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

# fit a multivariate BM with Penalized likelihood
fit <- fit_t_pl(Y, tree, model="BM", method="RidgeAlt")

# Perform the ancestral states reconstruction
ancestral <- estim(fit)

# retrieve the scores
head(ancestral$nodes)
}
}
