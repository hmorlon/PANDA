\name{plot_phylosignal_sub_network}
\alias{plot_phylosignal_sub_network}

\title{
Plot clade-specific phylogenetic signals in a bipartite interaction network
}

\description{
This function plots the clade-specific phylogenetic signals in species interactions. For each node of tree A having a certain number of descending species, it represents the phylogenetic signal in the resulting sub-network by performing a Mantel test between the phylogenetic distances and the ecological distances for the given sub-clade of tree A. 

}

\usage{
plot_phylosignal_sub_network(tree_A, results_sub_clades, network, legend=TRUE, show.tip.label=FALSE, where="bottomleft")
}

\arguments{

\item{tree_A}{
a phylogenetic tree of guild A (the columns of the interaction network). It must be an object of class "phylo".
}
\item{results_sub_clades}{
output of the function phylosignal_sub_network.
}

\item{network}{
a matrix representing the bipartite interaction network with species from guild A in columns and species from guild B in rows. Row names (resp. columns names) must correspond to the tip labels of tree B (resp. tree A).
}

\item{legend}{
indicates whether the legend should be plotted. 
}

\item{show.tip.label}{
indicates whether the tip labels should be plotted. 
}

\item{where}{
indicates where to put the legend (default is "bottomleft").
}
}

\value{
A phylogenetic tree with nodes colored according to the clade-specific phylogenetic signals. Blue nodes are not significant (Bonferonni correction), whereas orange-red nodes present significant phylogenetic signals and their color indicates the strength of the signal (correlation R of the Mantel test).
}
\details{

See the tutorial on GitHub (https://github.com/BPerezLamarque/Phylosignal_network). 

}
\references{
Perez‐Lamarque, B., Maliet, O., Selosse, M-A., Martos, F., Morlon H. (in prep). Do closely related species interact with similar partners? Testing for phylogenetic signal in bipartite interaction networks.

Goslee, S.C. & Urban, D.L. (2007). The ecodist package for dissimilarity-based analysis of ecological data. J. Stat. Softw., 22, 1–19.

Chen, J., Bittinger, K., Charlson, E.S., Hoffmann, C., Lewis, J., Wu, G.D., et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. Bioinformatics, 28, 2106–2113.



}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{phylosignal_network}}
\code{\link{phylosignal_sub_network}}
}
\examples{

# Load the data
data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # interaction matrix 
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)
tree_fungi <- mycorrhizal_network[[3]] # phylogenetic tree (phylo object)


# Clade-specific phylogenetic signal in species interactions in guild A 
# (do closely related species interact with similar partners in sub-clades of guild A?)

results_clade_A <- phylosignal_sub_network(network, tree_A = tree_orchids, tree_B = tree_fungi,
method = "GUniFrac", correlation = "Pearson")
plot_phylosignal_sub_network(tree_A = tree_orchids, results_clade_A, network)

# Clade-specific phylogenetic signal in species interactions in guild B 
# (do closely related species interact with similar partners in sub-clades of guild B?)

results_clade_B <- phylosignal_sub_network(t(network), tree_A = tree_fungi, tree_B = tree_orchids,
method = "GUniFrac", correlation = "Pearson")
plot_phylosignal_sub_network(tree_A = tree_fungi, results_clade_B, t(network))
}

