\name{fit_t_env}
\alias{fit_t_env}

\title{
Maximum likelihood fit of the environmental model of trait evolution
}
\description{
Fits model of trait evolution for which evolutionary rates depends on an environmental function.
}
\usage{

fit_t_env(phylo, data, env_data, error=NULL, model=c("EnvExp", "EnvLin"), par=NULL, 
       method="Nelder-Mead", control=list(maxit=20000))
}
\arguments{
  \item{phylo}{
An object of class 'phylo' (see ape documentation)
}
  \item{data}{
A named vector of phenotypic trait with the same size of the tree.
}
  \item{env_data}{
Environmental data, given as a time continuous function (see, e.g. splinefun) or a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).
}
  \item{error}{
A named vector with standard error for each species.
}
  \item{model}{
The model describing the functional form of variation of the evolutionary rate \eqn{\sigma} with time and the environmental variable. Default models are "EnvExp" and "EnvLin" (see details). An user defined function of any functional form may be used.
This function has three arguments: the first argument is time; the second argument is the environmental variable; the third arguement is a numeric vector of the parameters controlling the time and environmental variation (to be estimated). See the example below.
}
  \item{par}{
List of arguments to be passed to the function. See details. 
}
  \item{method}{
Methods used by the optimization routine (see ?optim for details).
}
  \item{control}{
Max. bound for the number of iteration of the optimizer; other options can be fixed on the list (see ?optim). 
}
}


\value{
	a list with the following components
  
  \item{LH}{the maximum log-likelihood value}
  \item{aic}{the Akaike's Information Criterion}
  \item{aicc}{the size sample corrected Akaike's Information Criterion}
  \item{free.parameters}{the number of estimated parameters}
  \item{param}{a numeric vector of estimated parameters, sigma and beta respectively for the defaults models. In the same order as defined by the user if a     customized model is provided}
  \item{root}{the estimated root value}
  \item{convergence}{convergence status of the optimizing function; "0" indicates convergence (See ?optim for details)}
  \item{hess.value}{reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached}
  \item{env_func}{the environmental function}
  \item{tot_time}{the root age of the tree}
  \item{model}{the fitted model (default models or user specified)}
}

\details{The various options are provided through the "par" list. 

-param: The starting values used for the model. The lengths of the "param" argument in the par list is used to compute the total number of parameters for user specified models.

-scale: scale the amplitude of the environmental curve between 0 and 1. This improve the parameters search.

-df: the degree of freedom to use to define the spline. As a default, smooth.spline(env_data[,1], env_data[,2])$df is used. See \emph{sm.spline} for details.

-upper: the upper bound for the parameter search when the "L-BFGS-B" method is used. See \emph{optim} for details.

-lower: the lower bound for the parameter search when the "L-BFGS-B" method is used. See \emph{optim} for details.

-sig2: can be used instead of param to define the starting sigma value only

-beta: can be used instead of param to define the beta starting value only
%The rate models here are function of environmental data fixed by
%env_data.
%Hence, if the parameter controlling the variation of \eqn{\lambda} with time is estimated to be positive (for example), this means that the speciation rate decreases from past to present.
}

\note{
	The speed of convergence of the fit might depend on the degree of freedom chosen to define the spline. 
}
\references{
Clavel, J. & Morlon, H. (in prep) 
}

\author{J. Clavel}

\seealso{
\code{\link{plot.fit_t.env}}, \code{\link{simulate_t_env}}, \code{\link{likelihood_t_env}}
}

\examples{
data(Cetacea)
data(InfTemp)

# Simulate trait on the Cetacean tree
set.seed(123)
trait<-rTraitCont(Cetacea)

## Fit the Environmental-exponential model
  # Fit the environmental model
  result1=fit_t_env(Cetacea, trait, env_data=InfTemp, par=list(scale=TRUE))
  plot(result1)

  # Add to the plot the results from different smoothing of the temperature curve
  result2=fit_t_env(Cetacea, trait, env_data=InfTemp, par=list(df=10, scale=TRUE))
  lines(result2, col="red")

  result3=fit_t_env(Cetacea, trait, env_data=InfTemp, par=list(df=50, scale=TRUE))
  lines(result3, col="blue")

## Fit the environmental linear model

  fit_t_env(Cetacea, trait, env_data=InfTemp, model="EnvLin", par=list(df=50, scale=TRUE))

## Fit user defined model

  # We define the function for the model
  my_fun<-function(t, env_cont, param){ 
      exp(param[1])*exp(param[2]*env_cont(t))
  }

  res<-fit_t_env(Cetacea, trait, env_data=InfTemp, model=my_fun, par=list(param=c(0.1,0)))
  res
  plot(res, col="red")

## Fit user defined environmental function

  spline_result <- sm.spline(InfTemp[,1],InfTemp[,2], df=50)
  env_func <- function(t){predict(spline_result,t)}
  t<-unique(InfTemp[,1])
  # We build the interpolated smoothing spline function
  env_data<-splinefun(t,env_func(t))
  
  # We then fit the model
  fit_t_env(Cetacea, trait, env_data=env_data)
  
## Various parameterization (box constraints, df, scaling of the curve...) example
# fit_t_env(Cetacea, trait, env_data=InfTemp, model="EnvLin", method="L-BFGS-B", par=list(scale=TRUE, lower=-30, upper=20, df=10))
}
