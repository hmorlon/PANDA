\name{createModelCoevolution}
\alias{createModelCoevolution}
\title{
Creation of a PhenotypicGMM
}
\description{
Creates an object of class PhenotypicGMM, a subclass of the class PhenotypicModel intended to represent the Generalist Matching Mutualism model of trait evolution on two specific trees.
}
\usage{
createModel(tree1, tree2)
}

\arguments{
\item{tree}{
an object of class 'phylo' as defined in the R package 'ape'.
}
\item{tree2}{
an object of class 'phylo' as defined in the R package 'ape'.
}
}

\value{
an object of class "PhenotypicGMM".
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

