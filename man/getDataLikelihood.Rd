\name{getDataLikelihood}
\alias{getDataLikelihood}
\title{
Likelihood of tip trait values.
}
\description{
Computes -log( likelihood ) of tip trait data under a given set of parameters, and for a specified model of trait evolution.
}
\usage{
getTipDistribution(model, data, params, v)
}

\arguments{
\item{model}{an object of class 'PhenotypicModel'.}
\item{data}{vector of tip trait data.}
\item{params}{vector of parameters, given in the same order as in the 'model' object.}
\item{v}{boolean specifying the verbose mode. Default value : FALSE.}
}

\value{
A numerical value : -log( likelihood )
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

