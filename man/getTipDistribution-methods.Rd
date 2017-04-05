\name{getTipDistribution-methods}
\docType{methods}
\alias{getTipDistribution-methods}
\alias{getTipDistribution,PhenotypicACDC-method}
\alias{getTipDistribution,PhenotypicADiag-method}
\alias{getTipDistribution,PhenotypicBM-method}
\alias{getTipDistribution,PhenotypicDD-method}
\alias{getTipDistribution,PhenotypicGMM-method}
\alias{getTipDistribution,PhenotypicModel-method}
\alias{getTipDistribution,PhenotypicOU-method}
\alias{getTipDistribution,PhenotypicPM-method}
\title{ Distribution of tip trait values. }
\description{
 Computes the mean and variance of the tip trait distribution under a specified model of trait evolution.
}

\usage{getTipDistribution(model, params, v)}

\arguments{
\item{model}{an object of class 'PhenotypicModel'}
\item{params}{vector of parameters, given in the same order as in the 'model' object.}
\item{v}{boolean specifying the verbose mode. Default value : FALSE.}
}

\value{
\item{mean}{Expectation vector of the tip trait distribution.}
\item{Sigma}{Variance-covariance matrix of the tip trait distribution.}
}

\section{Methods}{
\describe{

\item{\code{signature(object = "PhenotypicModel")}}{
In the most general case, this function computes the expectation vector and the variance-covariance matrix using a numerical integration procedure that may take time.
}

\item{\code{signature(object = "PhenotypicACDC")}}{
The function has been optimized for this subclass.
}

\item{\code{signature(object = "PhenotypicADiag")}}{
The function has been optimized for this subclass.
}

\item{\code{signature(object = "PhenotypicBM")}}{
The function has been optimized for this subclass.
}

\item{\code{signature(object = "PhenotypicDD")}}{
The function has been optimized for this subclass.
}

\item{\code{signature(object = "PhenotypicGMM")}}{
The function has been optimized for this subclass.
}

\item{\code{signature(object = "PhenotypicOU")}}{
The function has been optimized for this subclass.
}

\item{\code{signature(object = "PhenotypicPM")}}{
The function has been optimized for this subclass.
}
}}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
