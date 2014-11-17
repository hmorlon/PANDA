\name{print.fit.bd}
\alias{print.fit.bd}

\title{
Print fit_bd results
}
\description{
The function prints a compact summary of a fitted object result on the console.
}

\usage{
\method{print}{fit.bd}(x, \dots)
}

\arguments{
  \item{x}{
an object result of 'fit_bd'.
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

# Fit the pure birth model (no extinction) with exponential variation of the speciation rate with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.08, 0.01)
mu_par<-c()
result <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89,
                     expo.lamb = TRUE, fix.mu=TRUE)
# Generic print
print(result)
}
