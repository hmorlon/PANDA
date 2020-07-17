\name{delineate_phylotypes}
\alias{delineate_phylotypes}

\title{
Delineate phylotypes
}

\description{
This function traverses a tree from the root to the tips, at every node computes the average similarity of all sequences descending from the node, and collapses the sequences into a single phylotype if their sequence dissimilarity is lower than a given threshold. The average similarity can be computed using raw measured of the average similarity or using measures of genetic diversity (nucleotidic diveristy "pi" or Watterson "theta") which correct for gaps in the nucleotidic alignments (Ferretti et al., 2012). 
}

\usage{
delineate_phylotypes(tree, thresh=97, sequences, method="pi")
}

\arguments{

\item{tree}{
a phylogenetic tree of all the sequences. It must be an object of class "phylo" and must be rooted.
}
\item{thresh}{
a numeric digit between 0 and 100 indicating the minimal average similarity to collapse sequences within the same phylotype. By default the average similarity is 97.
}

\item{sequences}{
a matrix representing the nucleotidic alignment of all the sequences present in the phylogenetic tree. 
}

\item{method}{
indicates which method to use to compute the average similarity: "mean" computes the average raw distances between pairs of sequences, "pi" (default) measures the nucleotidic diversity (Nei & Li, 1979) while controlling for gaps in the alignment, and "theta" measures the Watterson theta genetic diveristy  (Watterson, 1975) also controlling for gaps. 
}

}


\value{
A table with row names corresponding to the sequence names. The first column to the phylotype assignation and the second colums indicates the name of the representative sequence of each phylotype (longuest sequence available). Phylotypes are numbered starting at 1, and all the phylotypes nammed "0" correspond to singletons.
}

\references{
Perez-Lamarque B, Öpik M, Maliet O, Silva A, Selosse MA, Martos F, and Morlon H (in prep). Global drivers of obligate mycorrhizal symbionts diversification

Ferretti L, Raineri E, Ramos-Onsins S. 2012. Neutrality tests for sequences with missing data. Genetics 191: 1397–1401.

Morlon H, O’Connor TK, Bryant JA, Charkoudian LK, Docherty KM, Jones E, Kembel SW, Green JL, Bohannan BJM. 2015. The biogeography of putative microbial antibiotic production. PLoS ONE 10.

Nei M. & Li W.-H., Mathematical model for studying genetic variation in terms of re- striction endonucleases, 1979, Proc. Natl. Acad. Sci. USA.

Watterson G.A. , On the number of segregating sites in genetical models without recombi-
nation, 1975, Theor. Popul. Biol.

}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{pi_estimator}}
\code{\link{theta_estimator}}
}
\examples{


}

