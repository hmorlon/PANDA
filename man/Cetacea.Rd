\name{Cetacea}
\alias{Cetacea}
\docType{data}
\title{Cetacean phylogeny}
\description{Ultrametric phylogenetic tree for 87 of the 89 extant cetacean species}
\usage{data(Cetacea)}
\format{
List with the following 4 components:
\itemize{
  \item edge, a two-column matrix of mode numeric where each row
    represents an edge of the tree; the nodes and the tips are
    symbolized with numbers; the tips are numbered 1, 2, \dots, and the
    nodes are numbered after the tips. For each row, the first column
    gives the ancestor.
  \item edge.length, a numeric vector giving the lengths of the
    branches given by \code{edge}.
  \item tip.label, a vector of mode character giving the names of the
    tips; the order of the names in this vector corresponds to the
    (positive) number in \code{edge}.
  \item Nnode, the number of (internal) nodes.}
% - attr(*, "class")= chr "phylo"
% - attr(*, "order")= chr "cladewise"
}
\details{
Phylogeny constructed by Steeman et al (see reference), which consists of 87 of 89 extant cetaceans species.
This phylogeny was derived from six mitochondrial and nine nuclear genes using the Bayesian phylogenetic inference.
It was calibrated using seven paleontological age constraints and a relaxed molecular clock approch.
}

\source{
Steeman ME et al.(2009) Radiation of extant cetaceans driven by restructuring of the oceans, Syst Biol 58:573-585
}

\references{
Steeman ME et al.(2009) Radiation of extant cetaceans driven by restructuring of the oceans, Syst Biol 58:573-585
Morlon, H., Parsons, T.L. and Plotkin, J. (2011) Reconciling molecular phylogenies with the fossil record, Proc Nat Acad Sci 108: 16327-16332
Condamine, F.L., Rolland, J., and Morlon, H. (2013) Macroevolutionary perspectives to environmental change, Eco Lett 16: 72-85
}
\examples{
data(Cetacea)
print(Cetacea)
plot(Cetacea)
}
\keyword{datasets}
