\name{likelihood_t_env}
\alias{likelihood_t_env}

\title{
Likelihood of a dataset under environmental models of trait evolution.
}
\description{
Computes the likelihood of a dataset under either the linear or exponential environmental model, or an user defined environmental model. This function is used internally by fit_t_env.
}

\usage{
likelihood_t_env(phylo, data, par=NULL, model=c("EnvExp", "EnvLin"))
}

\arguments{
  \item{phylo}{  
an object of class 'phylo' (see ape documentation)
}
  
  \item{data}{
a named vector of continuous data with names corresponding to \code{phylo$tip.label}
}

 \item{par}{
a list with "param", "fun", "times", "mtot" and "error" arguments.

-param: a vector with the parameters used in the environmental function. The first value is \code{sig2} and the second is \code{beta}.

-fun: a time contnuous function of an environmental variable (see e.g. ?fit_t_env)

-times: a vector of branching times starting at zero 
(e.g. max(branching.times(phylo))-branching.times(phylo))

-mtot: root age of the tree (e.g. max(branching.times(phylo)))

-error: a vector of standard error (se) for each species

If the "times" argument is not provided, the "phylo" object is used to compute it as well as "mtot".

Note that the argument "mu" can be used to specify the root state (e.g. when using an mcmc sampler)
}

 \item{model}{
model chosen to fit trait data, \code{"EnvExp"} is the exponential-environmental model, and \code{"EnvLin"} is the linear-environmental model. Otherwise, an user specified model can be provided.
}

}

\details{
	the "fun" argument can be filled by an environmental dataframe.
}

\value{
  the log-likelihood value of the environmental model
}
\references{
	Clavel, J., and Morlon, H. (in prep) 
}

\author{Julien Clavel
}

\seealso{
\code{\link{fit_t_env}}

}


\examples{
data(Cetacea)
data(InfTemp)

# Simulate trait on the Cetacean tree
set.seed(123)
trait<-rTraitCont(Cetacea)

# Compute the likelihood 
parameters <- list(param=c(0.1, 0), fun=InfTemp)
likelihood_t_env(Cetacea,trait,par=parameters,model="EnvExp")

# Provide the times
brtime<-branching.times(Cetacea)
mtot<-max(brtime)
times<-mtot-brtime

parameters <- list(param=c(0.1, 0), fun=InfTemp, times=times, mtot=mtot)
likelihood_t_env(Cetacea,trait,par=parameters,model="EnvExp")

# Provide the environmental function rather than the dataset (faster if used recursively)
\donttest{
spline_result <- sm.spline(InfTemp[,1],InfTemp[,2], df=50)
env_func <- function(t){predict(spline_result,t)}
t<-unique(InfTemp[,1])
# We build the interpolated smoothing spline function
env_data<-splinefun(t,env_func(t))
  
parameters <- list(param=c(0.1, 0), fun=env_data, times=times, mtot=mtot)
likelihood_t_env(Cetacea,trait,par=parameters,model="EnvExp")
	}  
}