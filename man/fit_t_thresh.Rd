\name{fit_t_thresh}
\alias{fit_t_thresh}

\title{
Maximum likelihood fit of Felsenstein' threshold model of evolution for a binary character.
}
\description{
Fits the threshold model of trait evolution for a binary (discrete) character. It uses various processes (Brownian motion, OU, EB, Pagel's lambda, and the ClimEnv model) for the latent variable that control the evolution of the binary character.
}
\usage{

fit_t_thresh(phylo, data, model=c("BM","OU","EB","lambda","Clim"), N=200,
          env_data=NULL, ...)
          
}
\arguments{
  \item{phylo}{
An object of class 'phylo' (see ape documentation)
}
  \item{data}{
A named vector of phenotypic trait values: a binary character (e.g. 0 and 1).
}

  \item{model}{
 Use a model for the continuous latent variable that controls the transitions between character states depending on the threshold value. Current models are "BM" (Brownian motion), "OU" (Ornstein Uhlenbeck), "EB" (Early Burst), "lambda" (Pagels' lambda), "Clim" (The ClimExp model implemented in the fit_t_env function - Clavel & Morlon 2017).
 }
 
  \item{N}{
N = is the number of knots used to approximate the integral (through numerical integration) in the calculation of the partial likelihoods. See Hiscott et al. (2016). Higher values should lead to more accurate approximations of the log-likelihood but will computationally more expensive.
}

  \item{env_data}{
Environmental data, given as a time continuous function (see, e.g. splinefun) or a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance). See also details in ?fit_t_env
}


\item{...}{
Optional arguments to be passed to the function.  
}
}


\value{
	a list with the following components
  
  \item{LH}{the maximum log-likelihood value}
  \item{aic}{the Akaike's Information Criterion}
  \item{aicc}{the second order Akaike’s Information Criterion}
  \item{free.parameters}{the number of estimated parameters}
  \item{param}{a numeric vector of estimated parameters for each models}
  \item{convergence}{convergence status of the optimizing function; "0" indicates convergence (See ?optim for details)}
  \item{env_func}{the environmental function - when the Clim model is used}
  \item{tot_time}{the root age of the tree - when the Clim model is used}
  \item{model}{the fitted model }
}

\details{

The \code{fit_t_thresh} function allows the threshold model (Felsenstein, 2012) to be fitted to a binary character by maximum likelihood using the pruning algorithm described by Hiscott et al. (2016). \code{fit_t_thresh} assumes that discrete trait evolution depends on an unseen continuous variable (referred to as 'liability') that evolves across the tree according to various processes ('Brownian motion', 'Ornstein-Uhlenbeck', 'Early-Burst', 'Pagel's lambda' and 'ClimEnv'; see Melendez-Vazquez et al. (2025) to test the effect of climate on the evolution of a binary character state (endothermy)). When using the 'Clim' model, this corresponds to the discrete trait counterpart of the 'ClimExp' model in fit_t_env (Clavel & Morlon, 2017).

The model estimates the parameters that control each evolutionary process for the latent traits (i.e. 'alpha' for OU, 'lambda' for Pagel's lambda, the deceleration rate 'a' for the EB model and the 'beta' parameter of the ClimExp model). It also estimates the value at the root for the latent process ('mu'). Note that, as the scale of the latent process ('sigma^2') and the threshold value are not identifiable separately, the value of the 'sigma^2' parameter in the above continuous processes is set to one.
}

\note{
	For the "Clim" model, the user defined function is evaluated forward in time i.e.: from the root to the tips (time = 0 at the (present) tips).
	The speed of convergence of the fit might depend on the degree of freedom chosen to define the spline. 
}

\references{
Felsenstein, J. 2012. A comparative method for both discrete and continuous characters using the threshold model. Am Nat. 179(2):145–156.

Hiscott, G., Fox, C., Parry, M., Bryant, D., 2016. Efficient recycled algorithms for quantitative trait models on phylogenies. Genome Biol. Evol. 8(5):1338-1350.

Clavel, J. & Morlon, H., 2017. Accelerated body size evolution during cold climatic periods in the Cenozoic. Proceedings of the National Academy of Science, 114(16): 4183-4188.

Melendez-Vazquez, F., Lucaci, A. G., Selberg A., Clavel, J., Rincon-Sandoval, M., Santaquiteria, A., White, W. T., Drabeck, D., Carnevale, G., Duarte-Ribeiro, E., Miya, M., Westneat, M. W., Baldwin, C. C., Hughes, L. C., Ortí, G., Kosakovsky Pond, S. L., Betancur-R, R., Arcila, D., 2025. Ecological interactions and genomic innovation fueled the evolution of ray-finned fish endothermy. Sci. Adv., 11 (eads8488):1-16.

}

\author{J. Clavel}

\seealso{
%% add later  \code{\link{sim_t_env}}  
\code{\link{plot.threshML}}
}

\examples{
\donttest{
## Comparison between models
## require(phytools)
set.seed(1)
tree <- pbtree(n=150, scale=67)
tree <- reorder(tree, "postorder")
trait <- rTraitCont(tree)
th <- sapply(trait, threshState, thresholds = setNames(c(0,Inf), letters[1:2]))

# Data were simulated with BM liabilities - Comparison of AIC
bm_fit <- fit_t_thresh(tree, th, model="BM")
lambda_fit <- fit_t_thresh(tree, th, model="lambda")
ou_fit <- fit_t_thresh(tree, th, model="OU")
eb_fit <- fit_t_thresh(tree, th, model="EB")

bm_fit$aic
lambda_fit$aic
ou_fit$aic
eb_fit$aic


# Fit the climatic model
data(InfTemp)
clim_fit <- fit_t_thresh(tree, th, model="Clim", env_data=InfTemp)


# Now simulate under the climatic model
# I use an (arbitrary) parameter value for the climatic model of -0.2
set.seed(1)
beta = -0.2
# The (latent) trait is then simulated using the function "sim_t_env" (see also ?fit_t_env)
x <- sim_t_env(tree, param=c(0.1, beta), env_data=InfTemp, model="EnvExp", root.value=0,
               step=0.001, plot=TRUE)

# Use the "threshState" function to simulate the discrete traits from the continuous one 'x'.
th <- sapply(x, threshState, thresholds = setNames(c(0,Inf), letters[1:2]))


clim_fit2 <- fit_t_thresh(tree, th, model="Clim", env_data=InfTemp)
bm_fit2 <- fit_t_thresh(tree, th, model="BM")

# compare the model fit
clim_fit2$aic
bm_fit2$aic
 }
}
