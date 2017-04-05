\name{createModel}
\alias{createModel}
\title{
Creation of a PhenotypicModel
}
\description{
Creates an object of class PhenotypicModel, intended to represent a model of trait evolution on a specific tree. DIstinct keywords correspond to different models, using one phylogenetic tree.
}
\usage{
createModel(tree, keyword)
}

\arguments{
\item{tree}{
an object of class 'phylo' as defined in the R package 'ape'.
}
\item{keyword}{
a string specifying the model. Available models include "BM", "BM_from0", "BM_from0_driftless", "OU", "OU_from0", "ACDC", "DD", "PM", "PM_OUless".
}
}

\value{
the object of class "PhenotypicModel".
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

