\name{Birds}
\alias{Birds}
\docType{data}
\title{Bird phylogeny}
\description{
Phylogenetic tree with 9993 tips and 9992 internal nodes.
}
\usage{data(Birds)}
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
}

\examples{
data(Birds)
print(Birds)
plot(Birds)
}
\keyword{datasets}
