\name{simulateTipData}
\alias{simulateTipData}
\title{
Tip trait simulation under a model of phenotypic evolution.
}
\description{
Simulates tip trait data under a specified model of phenotypic evolution, with three distinct behaviours specified with the 'method' argument.
}
\usage{
modelSelection(model, params, method, v)
}

\arguments{
\item{model}{an object of class 'PhenotypicModel'.}
\item{params}{vector of parameters, given in the same order as in the 'model' object.}
\item{method}{an integer specifying the behaviour of the function. If method = 1 (default value), the tip distribution is first computed, before returning a simulated dataset drawn in this distribution. If method = 2, the whole trajectory is simulated step by step, plotted, and returned. Otherwise, the whole trajectory is simulated step by step, and then returned without being plotted.}
}

\value{
a vector of trait values at the tips of the tree.
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

