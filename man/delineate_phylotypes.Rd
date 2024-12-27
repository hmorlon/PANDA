\name{delineate_phylotypes}
\alias{delineate_phylotypes}

\title{
Automatic phylotypes delineation 
}

\description{
This function traverses a tree from the root to the tips, at every node computes the average similarity of all sequences descending from the node, and collapses the sequences into a single phylotype if their sequence dissimilarity is lower than a given threshold. The average similarity can be computed using raw measured of the average similarity or using measures of genetic diversity (nucleotidic diversity "pi" (Nei & Li, 1979) or Watterson "theta" (Watterson, 1975)) which correct for gaps in the nucleotidic alignments (Ferretti et al., 2012). 
}

\usage{
delineate_phylotypes(tree, thresh=97, sequences, method="pi")
}

\arguments{

\item{tree}{
a phylogenetic tree of all the sequences. It must be an object of class "phylo" and must be rooted.
}
\item{thresh}{
a numeric digit between 0 and 100 indicating the minimal average similarity to collapse sequences within the same phylotype. By default, the average similarity is 97.
}

\item{sequences}{
a matrix representing the nucleotidic alignment of all the sequences present in the phylogenetic tree. 
}

\item{method}{
indicates which method to use to compute the average similarity: "mean" computes the average raw distances between pairs of sequences, "pi" (default) measures the nucleotidic diversity (Nei & Li, 1979) while controlling for gaps in the alignment, and "theta" measures the Watterson theta genetic diversity (Watterson, 1975) also controlling for gaps. 
}

}


\value{
A table with its row names corresponding to the sequence names. The first column corresponds to the phylotype assignation and the second columns indicates the name of the representative sequence of each phylotype (longest sequence available). Phylotypes are numbered starting at 1, and all the phylotypes named "0" correspond to singletons.
}

\references{
Perez-Lamarque B, Öpik M, Maliet O, Silva A, Selosse M-A, Martos F, and Morlon H. 2022. Analysing diversification dynamics using barcoding data: The case of an obligate mycorrhizal symbiont, Molecular Ecology, 31:3496–512.

Ferretti L, Raineri E, Ramos-Onsins S. 2012. Neutrality tests for sequences with missing data. Genetics 191: 1397–1401.

Morlon H, O’Connor TK, Bryant JA, Charkoudian LK, Docherty KM, Jones E, Kembel SW, Green JL, Bohannan BJM. 2015. The biogeography of putative microbial antibiotic production. PLoS ONE 10.

Nei M & Li WH, Mathematical model for studying genetic variation in terms of restriction endonucleases, 1979, Proc. Natl. Acad. Sci. USA.

Watterson GA , On the number of segregating sites in genetical models without recombination, 1975, Theor. Popul. Biol.

}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{pi_estimator}}
\code{\link{theta_estimator}}
}
\examples{

library(phytools)

data(woodmouse)

alignment <- as.character(woodmouse) # nucleotidic alignment 

tree <- midpoint.root(nj(dist.dna(woodmouse, pairwise.deletion = TRUE, 
model = "K80"))) # rooted neighbor-joining tree

\donttest{
delineate_phylotypes(tree, thresh = 99, alignment, method = "pi")
}

}

