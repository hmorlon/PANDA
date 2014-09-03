\name{plot_dtt}
\alias{plot_dtt}

\title{
Plot diversification function through time
}
\description{
Plot diversification function through time
}
\usage{
plot_dtt(fit.bd, tot_time, N0=1)
}

\arguments{
  \item{fit.bd}{
an object result of 'fit_bd'.
}
  \item{tot_time}{
the age of the underlying phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
  \item{N0}{
Number of extant species..
}

}


\value{
A list that countains:

  \item{time}{time values from 0 to tot_time on which we evaluate the diversification}
  \item{value}{Diversification values on previous time}
}

\seealso{
\code{\link[PANDA]{fit_bd}}
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
# plot ==> no need to use plot_dtt
# as plot is a generic function
plot_dtt(result, tot_time, N0=89)
}
