\name{phylosignal_sub_network}
\alias{phylosignal_sub_network}

\title{
Compute clade-specific phylogenetic signals in a bipartite interaction network
}

\description{
This function computes the clade-specific phylogenetic signals in species interactions. For each node of tree A having a certain number of descending species, it computes the phylogenetic signal in the resulting sub-network by performing a Mantel test between the phylogenetic distances and the ecological distances for the given sub-clade of tree A. 
Mantel tests can be computed using quantified or binary networks, with the Jaccard, Bray-Curtis, or UniFrac ecological distances.

}

\usage{
phylosignal_sub_network(network, tree_A, tree_B,
    method = "Jaccard_weighted", 
    correlation = "Pearson", 
    nperm = 1000, minimum=10, degree=F)
}

\arguments{

\item{network}{
a matrix representing the bipartite interaction network with species from guild A in columns and species from guild B in rows. Row names (resp. columns names) must correspond to the tip labels of tree B (resp. tree A).
}
\item{tree_A}{
a phylogenetic tree of guild A (the columns of the interaction network). It must be an object of class "phylo".
}
\item{tree_B}{
(optional) a phylogenetic tree of guild B (the rows of the interaction network). It must be an object of class "phylo".
}

\item{method}{
indicates which method is used to compute the phylogenetic signal in species interactions using Mantel tests. You can choose "Jaccard_weighted" (default) for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), "Bray-Curtis" for computing the Bray-Curtis dissimilarity, or "GUniFrac" for computing the weighted (or generalized) UniFrac distances ("UniFrac_unweighted" to not take into account the interaction abundances).
}

\item{correlation}{
indicates which correlation (R) must be used in the Mantel test, among Pearson (default), Spearman, and Kendall correlations.
}

\item{nperm}{
a number of permutations to evaluate the significance of the Mantel test. By default, it equals 10,000, but this can be very long for the Kendall correlation.
}

\item{minimum}{
indicates the minimal number of descending species for a node in tree A to compute its clade-specific phylogenetic signal. 
}

\item{degree}{
if degree=TRUE, Mantel tests testing for phylogenetic signal in the number of partners are additionally performed in each sub-clade.
}
}


\value{

For Mantel tests, the function outputs a table where each line corresponds to a tested clade and which contains at least 8 columns: the name of the node ("node"), the number of species in the sub-clade A ("nb_A"), the number of species in guild B associated with the sub-clade A ("nb_B"), the Mantel correlation for guild A ("mantel_cor"), its associated upper p-value ("pvalue_upper"), its associated lower p-value ("pvalue_lower"), and the associated Bonferroni corrected p-values ("pvalue_upper_corrected" and "pvalue_lower_corrected"). 

"mantel_cor" indicates the strength of the phylogenetic signal in the sub-clade A. The upper p-value "pvalue_upper" indicates the significance of the phylogenetic signal in the sub-clade A. The lower p-value "pvalue_lower" indicates the significance of the anti-phylogenetic signal in the sub-clade A. Both Bonferroni p-values are corrected using the number of tested nodes. 
For instance, if "pvalue_upper_corrected"<0.05 for a given node, there is a significant phylogenetic signal in the corresponding sub-clade of A.


If degree=TRUE, it also indicates in each sub-clade, the phylogenetic signal in the number of partners ("degree_mantel_cor") and its significance with or without the Bonferroni correction ("degree_pvalue_upper", "degree_pvalue_lower" and "degree_pvalue_upper_corrected", "degree_pvalue_lower_corrected")

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
\code{\link{phylosignal_sub_network}}
\code{\link{plot_phylosignal_sub_network}}
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
method = "GUniFrac", correlation = "Pearson", degree=T)
plot_phylosignal_sub_network(tree_A = tree_orchids, results_clade_A)

# Clade-specific phylogenetic signal in species interactions in guild B 
# (do closely related species interact with similar partners in sub-clades of guild B?)

results_clade_B <- phylosignal_sub_network(t(network), tree_A = tree_fungi, tree_B = tree_orchids, 
method = "GUniFrac", correlation = "Pearson", degree=T)
plot_phylosignal_sub_network(tree_A = tree_fungi, results_clade_B)

}

