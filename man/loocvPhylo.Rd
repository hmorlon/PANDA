\name{loocvPhylo}
\alias{loocvPhylo}

\title{
Phylogenetic Leave-One-Out Cross-Validated (LOOCV) penalized likelihood for high-dimensional models of trait evolution
}
\description{
Fits high-dimensional model of trait evolution on trees through penalized likelihood. A phylogenetic LOOCV is used to estimate model parameters}
\usage{

loocvPhylo(Y, tree, model=c("BM", "OU", "EB", "lambda"), method=c("RidgeAlt", "RidgeArch", "RidgeAltapprox", "LASSO", "LASSOapprox"), targM=c("null", "Variance", "unitVariance"), REML=TRUE, up=NULL, low=NULL, tol=NULL, starting=NULL, SE=NULL, scale.height=TRUE)
  
}
\arguments{
  \item{Y}{
A matrix of phenotypic traits values (the variables are represented as columns)
}
  \item{tree}{
An object of class 'phylo' (see ape documentation)
}
  \item{model}{
The evolutionary model, "BM" is Brownian Motion, "OU" is Ornstein-Uhlenbeck, "EB" is Early Burst, and "lambda" is Pagel's lambda transformation.
}
  \item{method}{
The penalty method. "RidgeArch": Archetype (linear) Ridge penalty, "RidgeAlt": Quadratic Ridge penalty, "LASSO": Least Absolute Selection and Shrinkage Operator. "RidgeAltapprox" and "LASSOapprox" are fast approximations of the LOOCV for the Ridge quadratic and LASSO penalties
}
  \item{targM}{
The target matrix used for the Ridge regularizations. "null" is a null target, "Variance" for a diagonal unequal variance target, "unitVariance" for an equal diagonal target. Only works with "RidgeArch","RidgeAlt", and "RidgeAltapprox" methods.
}

  \item{REML}{
Use REML (default) or ML for estimating the parameters.
}
  \item{up}{
Upper bound for the parameter search of the evolutionary model (optional).
}

 \item{low}{
Lower bound for the parameter search of the evolutionary model (optional).
}

  \item{tol}{
minimum value for the regularization parameter. Singularities can occur with a zero value in high-dimensional cases. (default is NULL)}


\item{starting}{
Starting values for the parameter search (optional).}

\item{SE}{
Standard errors associated with values in Y. If TRUE, SE will be estimated.}

\item{scale.height}{
Whether the tree should be scaled to unit length or not. (default is TRUE)}
}


\value{
	a list with the following components
  
  \item{loocv}{the (negative) cross-validated penalized likelihood}
  \item{model.par}{the evolutionary model parameter estimates}
  \item{gamma}{the regularization/tuning parameter of the penalized likelihood}
  \item{scaled_tree}{the phylogenetic tree with branch length stretched to the model estimated parameters}
  \item{model}{the evolutionary model}
  \item{method}{the penalization method}
  \item{p}{the number of traits}
  \item{n}{the number of species}
  \item{targM}{the target used for Ridge Penalization}
  \item{R}{a list with the estimated evolutionary covariance matrix and it's inverse}
  \item{REML}{logical indicating if the REML (TRUE) or ML (FALSE) method has been used}
  \item{Y}{the input dataset}
  \item{SE}{the estimated standard error}
}

\details{
\code{loocvPhylo} allows fitting various multivariate evolutionary models to high-dimensional datasets (where the number of variables p is larger than n). Models estimates are more accurate than maximum likelihood methods. Models fit can be compared using the GIC criterion (see ?GIC).
}

\note{
	The LASSO is computationally intensive. Please wait!
	The Ridge penalties with "null" or "unitVariance" targets are rotation invariants}
\references{
Clavel, J. Aristide, L., Morlon, H., 2017. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods. submitted
}

\author{J. Clavel}

\seealso{
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

# fit the model
loocvPhylo(Y, tree, model="BM", method="RidgeAlt")

# try on rotated axis (using PCA)
trans <- prcomp(Y, center=FALSE)
loocvPhylo(trans$x, tree, model="BM", method="RidgeAlt")
}
}
