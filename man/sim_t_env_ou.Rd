\name{sim_t_env_ou}
\alias{sim_t_env_ou}

\title{
Recursive simulation (root-to-tip) of the OU environmental model 
}
\description{
Simulates datasets for a given phylogeny under the OU environmental model (see ?fit_t_env_ou)
}
\usage{

sim_t_env_ou(phylo, param, env_data, model, step=0.01, 
              plot=FALSE, sigma, alpha, theta0, ...)


}
\arguments{
  \item{phylo}{
An object of class 'phylo' (see ape documentation)
}
  \item{param}{
A numeric vector of parameters for the user-defined climatic model. For the OU-environmental model, there is only one parameters (beta). If a model fit object of class 'fit_t_env.ou' is provided, the ML parameters are used to generate new datasets.
}
  \item{env_data}{
Environmental data, given as a time continuous function (see, e.g. splinefun) or a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).
}

  \item{model}{
The model describing the functional form of variation of the evolutionary trajectory of the optimum "theta(t)" with time and the environmental variable (see details for default model). An user defined function of any functional form can be used (forward in time).
This function has four arguments: the first argument is time; the second argument is the environmental variable; the third argument is a numeric vector of the parameters controlling the time and environmental variation (to be estimated), and the fourth is the theta_0 value. See the example below.
}

  \item{step}{
This argument describe the length of the segments to simulate across for the phylogeny. The smaller is the segment, the greater is the accuracy of the simulation at the expense of the computation time. 
}
  \item{plot}{
If TRUE, the simulated process is plotted.
}

  \item{sigma}{
The "sigma" parameter of the OU process.
}

  \item{alpha}{
The "alpha" parameter of the OU process.
}

  \item{theta0}{
The "theta" parameter at the root of the tree (t=0).
}

  \item{...}{
Arguments to be passed through. For instance, "col" for plot=TRUE. 
}
}


\value{
A named vector with simulated trait values for \eqn{n} species in the phylogeny
}

\details{
The users defined function is simulated forward in time i.e.: from the root to the tips. The speed of the simulations might depend on the value used for the "step" argument. 
It's possible to estimate the traits with the MLE from another fitted object (see the example below).
}

\references{
Clavel, J. & Morlon, H., 2017. Accelerated body size evolution during cold climatic periods in the Cenozoic. Proceedings of the National Academy of Sciences,  114(16): 4183-4188.

Troyer, E., Betancur-R, R., Hughes, L., Westneat, M., Carnevale, G., White W.T., Pogonoski, J.J., Tyler, J.C., Baldwin, C.C., Orti, G., Brinkworth, A., Clavel, J., Arcila, D., 2022. The impact of paleoclimatic changes on body size evolution in marine fishes. Proceedings of the National Academy of Sciences, 119 (29), e2122486119.

Goswami, A. & Clavel, J., 2024. Morphological evolution in a time of Phenomics. EcoEvoRxiv, https://doi.org/10.32942/X22G7Q

}

\author{J. Clavel}

\seealso{
%% add later  \code{\link{fit_t_env}}  
\code{\link{plot.fit_t.env}}, 
\code{\link{fit_t_env}},
\code{\link{fit_t_env_ou}},
\code{\link{plot.fit_t.env.ou}}
}

\examples{

\donttest{
\dontshow{test = FALSE}
if(test){

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

# define a model - here we replicate the default model used in fit_t_env_ou
my_model <- function(t, env, param, theta0) theta0 + param[1]*env(t)

# simulate the traits
trait <- sim_t_env_ou(tree, sigma=sqrt(sim_sigma2), alpha=sim_alpha,
                      theta0=sim_theta, param=beta, model=my_model,
                      env_data=InfTemp, step=0.01, scale=TRUE, plot=TRUE)

## Fit the Environmental model (default)

result_fit <- fit_t_env_ou(phylo = tree, data = trait, env_data =InfTemp,  
                          method = "Nelder-Mead", df=50, scale=TRUE)
plot(result_fit)


# We can also use the results from the previous fit to simulate a new dataset
trait2 <- sim_t_env_ou(tree, param=result_fit, step=0.001, plot=TRUE)

result_fit2 <- fit_t_env_ou(phylo = tree, data = trait2, env_data =InfTemp, 
                            method = "Nelder-Mead", df=50, scale=TRUE)
result_fit2
}
}

}
