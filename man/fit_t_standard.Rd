\name{fit_t_standard}
\alias{fit_t_standard}

\title{
Fits standard models of trait evolution incorporating known and nuisance measurement error
}
\description{
Fits Brownian motion (BM), Ornstein-Uhlenbeck (OU), or early burst (EB) models of trait evolution to a given dataset and phylogeny.
}
\usage{
fit_t_standard(phylo, data, model=c("BM","OU","EB"), error=NULL, two.regime=FALSE, 
		method="Nelder-Mead", echo=TRUE, ...)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation); if \code{two.regime=TRUE}, this must be a simmap object from \code{\link{make.simmap}} with two regimes
}
  \item{data}{
a named vector of trait values with names matching \code{phylo$tip.label}
}
  \item{model}{
model chosen to fit trait data, \code{"BM"} is the Brownian motion model, \code{"OU"} is the Ornstein-Uhlenbeck model, and \code{"EB"} is the early burst model.
}

  \item{error}{
A named vector with standard error (SE) for each species (with names matching \code{"phylo$tip.label"}). Default is NULL, if NA, then the SE is estimated from the data (a nuisance parameter for unknown errors). Note: When standard error are provided the nuisance parameter is also estimated.
}
  \item{two.regime}{
if \code{TRUE}, fits a two-regime model
}
  \item{method}{
optimization method from \code{link{optim}}
}

  \item{echo}{
prints information to console during fit
}

  \item{...}{
Optional arguments. e.g. "upper=xx", "lower=xx" to specify bounds on the parameter search. "fixedRoot=TRUE" to use an OU model where the root state is assumed fixed (instead of sampled from the stationary distribution)
}

}

\value{
	a list with the following elements:
  \item{LH}{maximum log-likelihood value}
  \item{aic}{Akaike Information Criterion value}
  \item{aicc}{AIC value corrected for small sample size}
  \item{free.parameters}{number of free parameters from the model}
  \item{sig2}{maximum-likelihood estimate of \code{sig2} parameter}
  \item{alpha}{maximum-likelihood estimate of \code{alpha} parameter of OU model (see Note)}
  \item{r}{maximum-likelihood estimate of the slope parameter of early burst model}
  \item{z0}{maximum-likelihood estimate of \code{z0}, the value at the root of the tree}
  \item{nuisance}{maximum-likelihood estimate of \code{nuisance}, the unknown, nuisance contribution to measurement error (see details)}
  \item{convergence}{convergence diagnostics from \code{\link{optim}} function (see optim documentation)}
}

\details{

Note: if including known measurement error, the model fit incorporates this known error and, in addition, estimates an unknown, nuisance contribution to measurement error. The current implementation does not differentiate between the two, so, for instance, it is not possible to estimate the nuisance measurement error without providing the known, intraspecific error values.
 
}


\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}

		Julien Clavel 
}

\seealso{
\code{\link{fit_t_comp}}
\code{\link{sim_t_tworegime}}
}

\examples{
\dontshow{test=FALSE}
\donttest{
if(test){
data(Cetacea_clades)
data<-sim_t_tworegime(Cetacea_clades,pars=c(sig2=0.01,r1=-0.01,r2=-0.02),
	root.value=0,Nsegments=1000,model="EB")
error<-rep(0.05,length(Cetacea_clades$tip.label))
names(error)<-Cetacea_clades$tip.label

#Fit single-regime models
BM1.fit<-fit_t_standard(Cetacea_clades,data,model="BM",error,two.regime=FALSE)
OU1.fit<-fit_t_standard(Cetacea_clades,data,model="OU",error,two.regime=FALSE)
EB1.fit<-fit_t_standard(Cetacea_clades,data,model="EB",error,two.regime=FALSE)

#Now fit models that incorporate biogeography, NOTE these models take longer to fit
BM2.fit<-fit_t_standard(Cetacea_clades,data,model="BM",error,two.regime=TRUE)
OU2.fit<-fit_t_standard(Cetacea_clades,data,model="OU",error,two.regime=TRUE)
EB2.fit<-fit_t_standard(Cetacea_clades,data,model="EB",error,two.regime=TRUE)
  }
}

}
