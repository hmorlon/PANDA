\name{mycorrhizal_network}
\alias{mycorrhizal_network}
\docType{data}
\title{Mycorrhizal network from La Réunion island}
\description{Mycorrhizal intercation network between orchids and mycorrhizal fungi from La Réunion island (Martos et al., 2012) along with the reconstructed phylogenetic trees of the orchids and the fungal OTUs.}
\usage{data(mycorrhizal_network)}
\details{
These phylogenies were constructed by maximum likelihood inference from four plastid genes for the orchids and one nuclear gene for the fungi. See Martos et al. (2012) for details.
}

\source{
Martos, F., Munoz, F., Pailler, T., Kottke, I., Gonneau, C. & Selosse, M.-A. (2012). The role of epiphytism in architecture and evolutionary constraint within mycorrhizal networks of tropical orchids. Mol. Ecol., 21, 5098–5109.
}

\references{

Martos, F., Munoz, F., Pailler, T., Kottke, I., Gonneau, C. & Selosse, M.-A. (2012). The role of epiphytism in architecture and evolutionary constraint within mycorrhizal networks of tropical orchids. Mol. Ecol., 21, 5098–5109.

Perez‐Lamarque, B., Maliet, O., Selosse, M-A., Martos, F., Morlon H. (in prep). Do closely related species interact with similar partners? Testing for phylogenetic signal in bipartite interaction networks.

}
\examples{

data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # interaction matrix 
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)
tree_fungi <- mycorrhizal_network[[3]] # phylogenetic tree (phylo object)

}
\keyword{datasets}
