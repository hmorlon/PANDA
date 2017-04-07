\name{createModelCoevolution}
\alias{createModelCoevolution}
\title{
Creation of a PhenotypicGMM
}
\description{
Creates an object of class PhenotypicGMM, a subclass of the class PhenotypicModel intended to represent the Generalist Matching Mutualism model of trait evolution on two specific trees.
}
\usage{
createModelCoevolution(tree1, tree2, keyword)
}

\arguments{
\item{tree1}{
an object of class 'phylo' as defined in the R package 'ape'.
}
\item{tree2}{
an object of class 'phylo' as defined in the R package 'ape'.
}
\item{keyword}{
a string object. Defaut value "GMM" returns an object of class PhenotypicGMM, which takes advantage of faster distribution computation. Otherwise, a "PhenotypicModel" is returned, and the computation of the tip distribution will take much longer.
}
}

\value{
an object of class "PhenotypicModel" or "PhenotypicGMM".
}

\references{
Manceau M., Lambert A., Morlon H. (2017) A unifying comparative phylogenetic framework including traits coevolving across interacting lineages Systematic Biology
}

\author{
M Manceau
}

