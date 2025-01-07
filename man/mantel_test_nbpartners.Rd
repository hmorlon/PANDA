\name{mantel_test_nbpartners}
\alias{mantel_test_nbpartners}

\title{
Compute Mantel test
}

\description{
This function tests for phylogenetic signal in species interactions in guild  A using a Mantel test that keep constant the number of partners per species.  
}

\usage{
mantel_test_nbpartners(network, tree_A, tree_B = NULL, method="Jaccard_binary",
nperm = 1000, correlation = "Pearson")
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
indicates which method is used to compute the phylogenetic signal in species interactions. If you want to perform a Mantel test between the phylogenetic distances and some ecological distances (do closely related species interact with similar partners?), you can choose "Jaccard_weighted" (default) for computing the ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), "Bray-Curtis" for computing the Bray-Curtis dissimilarity, or "GUniFrac" for computing the weighted (or generalized) UniFrac distances ("UniFrac_unweighted" to not take into account the interaction abundances).
}

\item{nperm}{
a number of permutations to evaluate the significance of the correlation. By default, it equals 1000.
}

\item{correlation}{
indicates which correlation (R) must be used among Pearson (default) and Spearman correlations.
}

\item{verbose}{
if TRUE, enables printing of messages.
}
}


\value{

  \item{mantelr }{Mantel correlation (R).}
  \item{pval1 }{one-tailed p-value (null hypothesis: R <= 0).}
  \item{pval2 }{one-tailed p-value (null hypothesis: R >= 0).}
  \item{pval3 }{two-tailed p-value (null hypothesis: R = 0).}

}

\references{
Perez-Lamarque B, Maliet O, Pichon B, Selosse M-A, Martos F, Morlon H. 2022. Do closely related species interact with similar partners? Testing for phylogenetic signal in bipartite interaction networks. Peer Community Journal, Volume 2, article no. e59. doi : 10.24072/pcjournal.179. https://peercommunityjournal.org/articles/10.24072/pcjournal.179/

Mantel, N. 1967. The detection of disease clustering and a generalized regression approach. Cancer Research 27:209-220.

}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{phylosignal_network}}

\code{\link{phylosignal_sub_network}}

\code{\link{mantel_test}}
}
\examples{

# Measuring phylogenetic signal in species interactions using a Mantel test 
# with permutations keeping constant the number of partners per species

library(RPANDA)

# Load the data
data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # bipartite interaction matrix 
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)

donttest({
mantel_test_nbpartners(network, tree_orchids, method="Jaccard_weighted", 
correlation="Pearson",  nperm = 1000)
})    
}

