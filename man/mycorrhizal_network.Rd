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

Martos, F., Munoz, F., Pailler, T., Kottke, I., Gonneau, C. & Selosse, M.-A. (2012). The role of epiphytism in architecture and evolutionary constraint within mycorrhizal networks of tropical orchids. Molecular Ecology, 21, 5098–5109.

Perez-Lamarque B, Maliet O, Pichon B, Selosse M-A, Martos F, Morlon H. 2022. Do closely related species interact with similar partners? Testing for phylogenetic signal in bipartite interaction networks. bioRxiv, 2021.08.30.458192, ver. 6 peer-reviewed and recommended by Peer Community in Evolutionary Biology. https://doi.org/10.1101/2021.08.30.458192

}
\examples{

data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # interaction matrix 
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)
tree_fungi <- mycorrhizal_network[[3]] # phylogenetic tree (phylo object)

}
\keyword{datasets}
