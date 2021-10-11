\name{mantel_test}
\alias{mantel_test}

\title{
Compute Mantel test
}

\description{
This function computes a Mantel test between two dissimilarity matrices. The available correlations are Pearson, Spearman, and Kendall. 
}

\usage{
mantel_test(formula = formula(data), data = sys.parent(), 
nperm = 1000, correlation = "Pearson")
}

\arguments{

\item{formula}{
formula y ~ x describing the test to be conducted where y and x  are distance matrices.
}
\item{data}{
an optional data frame containing the variables in the model as columns of dissimilarities. By default, the variables are taken from the current environment.
}

\item{correlation}{
indicates which correlation (R) must be used among Pearson (default), Spearman, and Kendall correlations.
}

\item{nperm}{
a number of permutations to evaluate the significance of the correlation. By default, it equals 1000, but this can be very long for the Kendall correlation.
}


}


\value{

  \item{mantelr }{Mantel correlation (R).}
  \item{pval1 }{one-tailed p-value (null hypothesis: R <= 0).}
  \item{pval2 }{one-tailed p-value (null hypothesis: R >= 0).}
  \item{pval3 }{two-tailed p-value (null hypothesis: R = 0).}

}
\details{

This function is adapted from the function mantel in the R-package ecodist (Goslee & Urban, 2007).

}
\references{
Perez‐Lamarque, B., Maliet, O., Selosse, M-A., Martos, F., Morlon H. (in prep). Do closely related species interact with similar partners? Testing for phylogenetic signal in bipartite interaction networks.

Goslee, S.C. & Urban, D.L. (2007). The ecodist package for dissimilarity-based analysis of ecological data. J. Stat. Softw., 22, 1–19.

Mantel, N. 1967. The detection of disease clustering and a generalized regression approach. Cancer Research 27:209-220.

}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{phylosignal_network}}

\code{\link{phylosignal_sub_network}}
}
\examples{

# Measuring phylogenetic signal in species interactions using a Mantel test 
# (do closely related species interact with similar partners?)

library(RPANDA)

# Load the data
data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # bipartite interaction matrix 
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)

network <- network[,tree_orchids$tip.label]

ecological_distances <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
    
phylogenetic_distances <- cophenetic.phylo(tree_orchids)

mantel_test(as.dist(ecological_distances) ~ as.dist(phylogenetic_distances), 
correlation="Pearson",  nperm = 10000)
      
}

