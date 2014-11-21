\name{sim_sgd}
\alias{sim_sgd}
\title{
Simulation method of a phylogenetic tree under SGD model
}
\description{
Simulate a SGD arising from the SGD model with exponentially increasing metapopulation size. Notations follow Manceau et al. (submitted).
}
\usage{
sim_sgd(tau, b, d, nu)
}

\arguments{
\item{tau}{
the total time of the process, i.e. the height of the phylogeny.
}
\item{b}{
the (constant) birth rate of individuals in the model.
}
\item{d}{
the (constant) death rate of individuals in the model.
}
\item{nu}{
the (constant) mutation rate of individuals in the model.
}
}

\value{
a phylogenetic tree of class "phylo". See ape for documentation.
}

\references{
Manceau M., Lambert A., Morlon H. , submitted.
}

\examples{
b <- 1e6
d <- b-0.5
nu <- 0.6
tree <- sim_bd(10,b,d,nu)
plot(tree)
}
