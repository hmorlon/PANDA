\name{lines.fit_t.env.ou}
\alias{lines.fit_t.env.ou}

\title{
Add to a plot line segments joining the phenotypic evolutionary optimum through time estimated by the \emph{fit_t_env_ou} function
}
\description{
Plot estimated optimum as a function of the environmental data and time.
}
\usage{

\method{lines}{fit_t.env.ou}(x, steps = 100, ...)

}
\arguments{
  \item{x}{
an object of class 'fit_t.env.ou' obtained from a \emph{fit_t_env_ou} fit.
}
  \item{steps}{
the number of steps from the root to the present used to compute the optimum \eqn{\theta(t)} through time.
}
  \item{...}{
further arguments to be passed to \code{plot}. See ?\code{plot}.
}

}


\value{
	\code{lines.fit_t.env.ou} returns invisibly a list with the following components used to add the line segments to the current plot:
  
  \item{time_steps}{the times steps where the climatic function was evaluated to compute the reconstructed optimum. The number of steps is controlled through the argument \code{steps}.}
  \item{values}{the estimated optimum through time estimated at each \code{time_steps}}
  
}


\note{
	All the graphical parameters (see \code{\link{par}}) can be passed through (e.g. line type: \code{lty}, line width: \code{lwd}, color: \code{col} ...)
}
\references{
Clavel, J. & Morlon, H., 2017. Accelerated body size evolution during cold climatic periods in the Cenozoic. Proceedings of the National Academy of Sciences,  114(16): 4183-4188.

Troyer, E., Betancur-R, R., Hughes, L., Westneat, M., Carnevale, G., White W.T., Pogonoski, J.J., Tyler, J.C., Baldwin, C.C., Orti, G., Brinkworth, A., Clavel, J., Arcila, D., 2022. The impact of paleoclimatic changes on body size evolution in marine fishes. Proceedings of the National Academy of Sciences, 119 (29), e2122486119.

Goswami, A. & Clavel, J., 2024. Morphological evolution in a time of Phenomics. EcoEvoRxiv, https://doi.org/10.32942/X22G7Q
}

\author{J. Clavel}

\seealso{
%% add later  \code{\link{simulate_t_env}}  
\code{\link{plot.fit_t.env.ou}}, \code{\link{fit_t_env_ou}}
}

\examples{
\dontshow{test = FALSE}
if(test){
\donttest{
data(InfTemp)
set.seed(9999) # for reproducibility

# Let's start by simulating a trait under a climatic OU
beta = 0.6           # relationship to the climate curve
sim_theta = 4        # value of the optimum if the relationship to the climate 
# curve is 0 (this corresponds to an 'intercept' in the linear relationship used below)
sim_sigma2 = 0.025   # variance of the scatter = sigma^2
sim_alpha = 0.36     # alpha value = strength of the OU; quite high here...
delta = 0.001        # time step used for the forward simulations => here its 1000y steps
tree <- pbtree(n=200, d=0.3) # simulate a bd tree with some extinct lineages
root_age = 60        # height of the root (almost all the Cenozoic here)
tree$edge.length <- root_age*tree$edge.length/max(nodeHeights(tree)) 
# here - for this contrived example - I scale the tree so that the root is at 60 Ma

trait <- sim_t_env_ou(tree, sigma=sqrt(sim_sigma2), alpha=sim_alpha, theta0=sim_theta, param=beta, 
              env_data=InfTemp, step=0.01, scale=TRUE, plot=FALSE)

## Fit the Environmental model (default)

result1 <- fit_t_env_ou(phylo = tree, data = trait, env_data =InfTemp,
                        method = "Nelder-Mead", df=50, scale=TRUE)
plot(result1, lty=2)

result2 <- fit_t_env_ou(phylo = tree, data = trait, env_data =InfTemp, 
                        method = "Nelder-Mead", df=10, scale=TRUE)
lines(result2, col="red")
}
}

}
