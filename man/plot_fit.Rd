\name{plot_fit}
\alias{plot_fit}

\title{
Plot speciation, extinction & rate functions of a fitted model
}
\description{
Plot speciation, extinction & diversification functions of a fitted model
}
\usage{
plot_fit(fit.bd, tot_time)
}

\arguments{
  \item{fit.bd}{
an object result of 'fit_bd'.
}
  \item{tot_time}{
the age of the underlying phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
}


\value{
  One or several graphical windows 
on which the speciation, extinction & diversification functions are represented depending on the fitted model (extinction/no extinction)
}


\details{
The fitted object countains the informations that are quickly usuable to plot the fitted functions.
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
# plot fitted rates
#plot_fit(result, tot_time)
}
