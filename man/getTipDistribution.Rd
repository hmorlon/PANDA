\name{getTipDistribution}
\alias{getTipDistribution}
\title{
Distribution of tip trait values.
}
\description{
Computes the mean and variance of the tip trait distribution under a specified model of trait evolution.
}
\usage{
getTipDistribution(object, params, v)
}

\arguments{
\item{object}{
an object of class 'PhenotypicModel'
}
\item{params}{
vector of parameters, given in the same order as in the 'model' object.
}
\item{v}{
boolean specifying the verbose mode. Default value : FALSE.
}
}

\value{
\item{mean}{Expectation vector of the tip trait distribution.}
\item{Sigma}{Variance-covariance matrix of the tip trait distribution.}
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

