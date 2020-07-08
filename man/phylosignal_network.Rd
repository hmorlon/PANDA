\name{phylosignal_network}
\alias{phylosignal_network}

\title{
Compute phylogenetic signal in a bipartite interaction network
}

\description{
This function computes the phylogenetic signal in a bipartite interaction network, either the phylogenetic signal in species interactions (do closely related species interact with similar partners?) using Mantel tests or the phylogenetic signal in the degree of generalism (do closely related species interact with the same number of partners?) using Mantel tests or the Phylogenetic bipartite linear model (PBLM) from Ives and Godfray (2006). 
Mantel tests measuring phylogenetic signal in species interactions can be computed using quantified or binary networks, with the Jaccard or UniFrac ecological distances.

}

\usage{
phylosignal_network(network, tree_A, tree_B=NULL, 
method = "Jaccard_weighted", nperm = 10000, 
correlation = "Pearson", only_A = FALSE)
}

\arguments{

\item{network}{
a matrix representing the ecological interaction network with species from guild A in columns and species from guild B in rows. Row names (resp. columns names) must correspond to the tip labels of tree B (resp. tree A).
}
\item{tree_A}{
a phylogenetic tree of the guild A (the columns of the interaction network). It must be an object of class "phylo".
}
\item{tree_B}{
(optional) a phylogenetic tree of the guild B (the rows of the interaction network). It must be an object of class "phylo".
}

\item{method}{
indicates which method is used to compute the phylogenetic signal in species interactions. If you want to perform a Mantel test between the phylogenetic distances and some ecological distances, you can choose "Jaccard_weighted" (default) for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the generalized UniFrac distances ("UniFrac_unweighted" to not take into account the interaction abundances).

Conversely, if you want to evaluate the phylogenetic signal in the degree of generalism, you can choose "degree". Alternatively (not recommended), you can use the Phylogenetic Bipartite Linear Model "PBLM" (see Ives and Godfray, 2006) or "PBLM_binary" to not consider the abundances of the interactions.
}

\item{correlation}{
(optional) indicates which correlation must be used in the Mantel test, among Pearson (default), Spearman, and Kendall correlations. It only applies for the methods "Jaccard_weighted", "Jaccard_binary", "GUniFrac", "UniFrac_unweighted", or "degree".
}

\item{nperm}{
(optional) number of permutations to evaluate the significance of the Mantel test. By default, it equals 10 000, but this can be very long for the Kendall correlation. It only applies for the methods "Jaccard_weighted", "Jaccard_binary", "GUniFrac", "UniFrac_unweighted", or "degree".
}

\item{only_A}{
(optional) indicates whether the signal should be only computed for guild A (and not for guild B). By default, it is computed for both guilds if "tree_B" is provided. 
}


}


\value{

For Mantel tests, the function outputs a vector of up to 8 values: the number of species in guild A ("nb_A"), the number of species in guild B ("nb_B"), the Mantel correlation for guild A ("mantel_cor_A"), its associated upper p-value ("pvalue_high_A"), its associated lower p-value ("pvalue_low_A"), and (optional) the Mantel correlation for guild B ("mantel_cor_B"), its associated upper p-value ("pvalue_high_B"), abd its associated lower p-value ("pvalue_low_B"),

"mantel_cor_A" (or "mantel_cor_B") indicates the strength of phylogenetic signal in guild A (or B). The upper p-value "pvalue_high_A" (or "pvalue_high_B") indicates the significance of phylogenetic signal in guild A (or B). The lower p-value "pvalue_low_A" (or "pvalue_low_B") indicates the significance of anti-phylogenetic signal in guild A (or B).  

For the PBLM approach (Ives and Godfray, 2006), the function outputs a vector of 8 values: the number of species in guild A ("nb_A"), the number of species in guild B ("nb_B"), the phylogenetic signals in guilds A ("dA") and B ("dB"), the covariance of interaction matrix ("MSETotal"), the mean square error of the complete model ("MSEFull"), the mean square error of model run on star phylogenies ("MSEStar"), and the mean square error of the model assuming strict Brownian motion evolutions ("MSEBase"). The significance of the phylogenetic signal can be evaluated by comparing "MSEFull" and "MSEStar".

}
\details{

See tutorial on GitHub (XXX). 

}
\references{
Perez‐Lamarque, B., Maliet, O., Selosse, M-A., Martos, F., Morlon H. (in prep). Measuring phylogenetic signal in bipartite interaction networks: Do closely related species interact with similar partners?

Goslee, S.C. & Urban, D.L. (2007). The ecodist package for dissimilarity-based analysis of ecological data. J. Stat. Softw., 22, 1–19.

Ives, A.R. & Godfray, H.C.J. (2006). Phylogenetic analysis of trophic associations. Am. Nat., 168, E1–E14.

Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D., et al. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26, 1463–1464.

Chen, J., Bittinger, K., Charlson, E.S., Hoffmann, C., Lewis, J., Wu, G.D., et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. Bioinformatics, 28, 2106–2113.



}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{phylosignal_sub_network}}
}
\examples{

# Load the data
data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # interaction matrix 
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)
tree_fungi <- mycorrhizal_network[[3]] # phylogenetic tree (phylo object)


# Using Mantel tests: 

# Phylogenetic signal in the degree of generalism (do closely related species interact with the same number of partners?)

phylosignal_network(network, tree_A = tree_orchids, method = "degree", correlation = "Pearson", nperm = 10000) # for guild A
phylosignal_network(t(network), tree_A = tree_fungi, method = "degree", correlation = "Pearson", nperm = 10000) # for guild B


# Phylogenetic signal in species interactions (do closely related species interact with similar partners?)

phylosignal_network(network, tree_A = tree_orchids, tree_B = tree_fungi, method = "GUniFrac", correlation = "Pearson", nperm = 10000) # measured for both guilds


# Using PBLM (not recommended)

# phylosignal_network(network, tree_A = tree_orchids, tree_B = tree_fungi, method = "PBLM") # very slow 

}

