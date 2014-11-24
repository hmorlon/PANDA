\name{fit_sgd}
\alias{fit_sgd}
\title{
Maximum likelihood fit of the SGD model
}

\description{
Fits the SGD model with exponential growth of the metacommunity, by maximum likelihood. Notations follow Manceau et al. (submitted)
}

\usage{
fit_sgd(phylo, tot_time, par_init, f=1)
}

\arguments{
\item{phylo}{
an object of type 'phylo' (see ape documentation)
}
\item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
\item{par_init}{
a numeric vector of initial values for the parameters (b,d,nu) to be estimated (these values are used by the optimization algorithm). 
}
\item{f}{
the fraction of extant species included in the phylogeny
}
}

\value{
the vector with estimated values of b (birth), b-d (growth) and nu (mutation).
}
\references{
Manceau M., Lambert A., Morlon H. , submitted.
}

\seealso{
\code{\link[PANDA]{likelihood_sgd}}
}
\examples{
# Some examples may take a little bit of time. Be patient!
data(Cetacea)
tot_time <- max(node.age(Cetacea)$ages)
par_init <- c(1., 0.9, 0.5)
result <- fit_sgd(Cetacea, par_init, f=87/89)$par
}
