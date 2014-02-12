\name{Phocoenidae}
\alias{Phocoenidae}
\docType{data}
\title{Phocoenidae phylogeny}
\description{
Ultrametric phylogenetic tree of the 6 extant Phocoenidae (porpoise) species
}
\usage{data(Phocoenidae)}
%\format{
%List with the following 5 components:
%\itemize{
%  \item edge, a two-column matrix of mode numeric where each row
%    represents an edge of the tree; the nodes and the tips are
%    symbolized with numbers; the tips are numbered 1, 2, \dots, and the
%    nodes are numbered after the tips. For each row, the first column
%    gives the ancestor.
%  \item edge.length, a numeric vector giving the lengths of the
%    branches given by \code{edge}.
%  \item tip.label, a vector of mode character giving the names of the
%    tips; the order of the names in this vector corresponds to the
%    (positive) number in \code{edge}.
%  \item Nnode, the number of (internal) nodes.
%  \item root.edge: a numeric value giving the length of the branch
%          at the root if it exists.
%     If several trees are read in the file, the returned object is of
%     class "multiPhylo", and is a list of objects of class '"phylo"'.
%     The name of each tree can be specified by 'tree.names', or can be
%     read from the file (see details).
%}
%}
\details{This phylogeny was extracted from Steeman et al. Syst Bio 2009 cetacean phylogeny}

\references{
Steeman ME et al.(2009) Radiation of extant cetaceans driven by restructuring of the oceans \emph{Syst Biol} 58:573-585

Morlon, H., Parsons, T.L., Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
}

\examples{
data(Phocoenidae)
print(Phocoenidae)
plot(Phocoenidae)
}
\keyword{datasets}
