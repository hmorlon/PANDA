\name{print.fit.env}
\alias{print.fit.env}

\title{
Print fit_env results
}
\description{
The function prints a compact summary of a fitted object result on the console.
}

\usage{
\method{print}{fit.env}(x, \dots)
}

\arguments{
  \item{x}{
an object result of 'fit_env'.
}
  \item{...}{
further arguments passed to or from other methods.
}
}


\value{
NULL
}

\seealso{
\code{\link[base]{print}}
}

\examples{

data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
data(InfTemp)
dof<-smooth.spline(InfTemp[,1], InfTemp[,2])$df

# fits a model with lambda varying as an exponential function of temperature and mu fixed to 0 (no extinction)
# here t stands for time and x for temperature

f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10, 0.01)
mu_par<-c()
result <- fit_env(Cetacea,InfTemp,tot_time,f.lamb,f.mu,
      lamb_par,mu_par,f=87/89, fix.mu=TRUE, df=dof, dt=1e-3)

# Generic print
print(result)
}
