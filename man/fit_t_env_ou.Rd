\name{fit_t_env_ou}
\alias{fit_t_env_ou}

\title{
Maximum likelihood fit of the OU environmental model of trait evolution
}
\description{
Fits Ornstein-Uhlenbeck (OU) model of trait evolution for which the optimum depends on an environmental function, or more generally a time varying function.
}
\usage{

fit_t_env_ou(phylo, data, env_data, error=NULL, model,
          method="Nelder-Mead", control=list(maxit=20000), ...)
          
}
\arguments{
  \item{phylo}{
An object of class 'phylo' (see ape documentation)
}
  \item{data}{
A named vector of phenotypic trait values.
}
  \item{env_data}{
Environmental data, given as a time continuous function (see, e.g. splinefun) or a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).
}
  \item{error}{
A named vector with standard errors (SE) of trait values for each species (with names matching \code{"phylo$tip.label"}). The default is NULL, in this case potential error is ignored in the fit. If set to NA, the SE is estimated from the data (to be used when there are no error measurements, a nuisance parameter is estimated). Note: When standard errors are provided, a nuisance parameter is also estimated.
}

 \item{model}{
 A user defined model. If not provided, a default model is used (see details)
 }

  \item{method}{
Methods used by the optimization routine (see ?optim for details).
}
  \item{control}{
Max. bound for the number of iteration of the optimizer; other options can be fixed on the list (see ?optim). 
}

\item{...}{
Arguments to be passed to the function. See details. 
}
}


\value{
	a list with the following components
  
  \item{LH}{the maximum log-likelihood value}
  \item{aic}{the Akaike's Information Criterion}
  \item{aicc}{the second order Akaike’s Information Criterion}
  \item{free.parameters}{the number of estimated parameters}
  \item{param}{a numeric vector of estimated parameters, sigma and beta respectively for the defaults models. In the same order as defined by the user if a custom model is provided}
  \item{root}{the estimated root value}
  \item{convergence}{convergence status of the optimizing function; "0" indicates convergence (See ?optim for details)}
  \item{hess.value}{reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached}
  \item{env_func}{the environmental function}
  \item{tot_time}{the root age of the tree}
  \item{model}{the fitted model (default models or user specified)}
  \item{nuisance}{the estimated SE for species mean when "error=NA"}
}

\details{
\code{fit_t_env_ou} allows fitting OU-environmental models of trait evolution (Troyer et al. 2020, Goswami & Clavel 2024). Compared to model implemented in \code{fit_t_env} where the rate of phenotypic evolution evolves as a function of an environmental variable (Clavel & Morlon 2020), here it's the optimum of a generalized Ornstein-Uhlenbeck (also called Hull-White model) that can changes as a function of an environmental variable T(t). More formally, the model is defined by the following process:

\deqn{dX(t)  =  \alpha  (\theta(t) -X(t))dt + \sigma dB(t)}{dX(t) = alpha*(theta(t) - X(t))dt + sigma*dB(t) }

Note that this model works only on NON-ULTRAMETRIC trees (e.g., with fossils)


The default model has the optimum changing as a function of environmental changes though times as 
defined below:

\deqn{\theta (t)  =  \theta_0 + \beta T(t)}{theta(t) = theta_0 + beta * T(t)}

Users defined models should have the following form (see also examples below):

\code{fun <- function(t, env, param, theta0){ theta0 + param*env(t)}}

\emph{t}: is the time parameter.

\emph{env}: is a time function of an environmental variable. 
See for instance object created by \code{splinefun} when interpolating coordinate of points.

\emph{param}: is a vector of parameters to estimate.

\emph{theta_0}: is the state at the root of the tree.

For instance, the default model function can be coded as: 

\code{fun <- function(t, env, param, theta0){ theta0 + param[1]*env(t)}} 

where \code{param[1]} is the \eqn{\beta} parameter. 
Note that in this case, one starting value should be provided in the \code{param} 
argument.

e.g.:

\code{beta=0}

\code{fit_t_env(tree, data, env_data=InfTemp, model=fun, param=beta)}



The various options are passed through "...". 

-param: The starting values used for the model. Must match the total number of parameters of the specified models. If "error=NA", a starting value for the SE to be estimated must be provided with user-defined models.

-scale: scale the amplitude of the environmental curve between 0 and 1. This may improve the parameters search in some situations.

-df: the degree of freedom to use for defining the spline. As a default, smooth.spline(env_data[,1], env_data[,2])$df is used. See \emph{sm.spline} for details.

-upper: the upper bound for the parameter search when the "L-BFGS-B" method is used. See \emph{optim} for details.

-lower: the lower bound for the parameter search when the "L-BFGS-B" method is used. See \emph{optim} for details.

-maxdiff: difference in time between tips and present day for phylogenetic trees with no contemporaneous species (default is 0)
}

\note{
	The users defined function is evaluated forward in time i.e.: from the root to the tips (time = 0 at the (present) tips).
	The speed of convergence of the fit might depend on the degree of freedom chosen to define the spline. 
}
\references{

Clavel, J. & Morlon, H., 2017. Accelerated body size evolution during cold climatic periods in the Cenozoic. Proceedings of the National Academy of Science, 114(16): 4183-4188.

Troyer, E., Betancur-R, R., Hughes, L., Westneat, M., Carnevale, G., White W.T., Pogonoski, J.J., Tyler, J.C., Baldwin, C.C., Orti, G., Brinkworth, A., Clavel, J., Arcila, D., 2022 - The impact of paleoclimatic changes on body size evolution in marine fishes. Proceedings of the National Academy of Sciences, 119 (29), e2122486119.

Goswami, A. & Clavel, J., 2024. Morphological evolution in a time of Phenomics. EcoEvoRxiv, https://doi.org/10.32942/X22G7Q

}

\author{J. Clavel}

\seealso{
%% add later  \code{\link{sim_t_env}}  
\code{\link{plot.fit_t.env.ou}},\code{\link{sim_t_env_ou}}
}

\examples{

data(InfTemp)

# Simulate a trait with temperature dependence of the optimum on a simulated tree
\donttest{

set.seed(9999) # for reproducibility

# Let's start by simulating a trait under a climatic OU
beta = 0.6           # relationship to the climate curve
sim_theta = 4        # value of the optimum if the relationship to the climate curve is 0 
sim_sigma2 = 0.025   # variance of the scatter = sigma^2
sim_alpha = 0.36     # alpha value = strength of the OU; quite high here...
delta = 0.001        # time step used for the forward simulations => here its 1000y steps
tree <- phytools::pbtree(n=200, d=0.3) # simulate a bd tree with some extinct lineages
root_age = 60        # height of the root (almost all the Cenozoic here)
tree$edge.length <- root_age*tree$edge.length/max(phytools::nodeHeights(tree)) 
# here - for this contrived example - I scale the tree so that the root is at 60 Ma

trait <- sim_t_env_ou(tree, sigma=sqrt(sim_sigma2), alpha=sim_alpha, theta0=sim_theta, 
                      param=beta, env_data=InfTemp, step=0.01, scale=TRUE, plot=TRUE)

## Fit the Environmental model (default)

result1 <- fit_t_env_ou(phylo = tree, data = trait, env_data =InfTemp, 
                        method = "Nelder-Mead", df=50, scale=TRUE)
plot(result1)


## Fit user defined model (note that several other environmental variables 
## can be simultaneously encapsulated in this function through the env argument)

# We re-define the function for the OU model with linear trend to the climatic curve
# NOTE: the env(t) function should return the value at the root for t=0

my_fun<-function(t, env, param, theta0){ 
    theta0 + param[1]*env(t)
}
  
# starting value for param[1]. Here we use an arbitrary value of 0.1
beta_guess = 0.1

# fit the model
result2 <- fit_t_env_ou(phylo = tree, data = trait, env_data =InfTemp, 
                        model = my_fun, param = beta_guess,  
                        method = "Nelder-Mead", df=50, scale=TRUE)
                  
# Retrieve the parameters and compare to 'result1'
result2
lines(result2, col="red", lty=2)


## Fit user defined environmental function

require(pspline)
  	 spline_result <- sm.spline(x=InfTemp[,1],y=InfTemp[,2], df=50)
  	 env_func <- function(t){predict(spline_result,t)}
  	 t<-unique(InfTemp[,1])
  	
  # We build the interpolated smoothing spline function (not scaled here)
  	 env_data<-splinefun(t,env_func(t))
  
  # We then fit the model
  
result3 <- fit_t_env_ou(phylo = tree, data = trait, env_data = env_data, 
                        model = my_fun, param = 0.01, method = "Nelder-Mead")

 }
}
