\name{theta_estimator}
\alias{theta_estimator}

\title{
Compute Watterson genetic diversity (Theta estimator)
}

\description{
This function computes the Theta estimator of genetic diversity (Watterson, 1975) while controlling for the presence of gaps in the alignment (Ferretti et al, 2012), frequent in barcoding datasets. 
}

\usage{
theta_estimator(sequences)
}

\arguments{

\item{sequences}{
a matrix representing the nucleotidic alignment of all the sequences present in the phylogenetic tree. 
}

}


\value{
An estimate of genetic diversity. 
}

\references{
Watterson GA , On the number of segregating sites in genetical models without recombination, 1975, Theor. Popul. Biol.

Ferretti L, Raineri E, Ramos-Onsins S. 2012. Neutrality tests for sequences with missing data. Genetics 191: 1397–1401.

Perez-Lamarque B, Öpik M, Maliet O, Silva A, Selosse M-A, Martos F, and Morlon H (in prep). Global drivers of obligate mycorrhizal symbionts diversification
}
\author{
Ana C. Afonso Silva & Benoît Perez-Lamarque}
\seealso{
\code{\link{pi_estimator}}
\code{\link{delineate_phylotypes}}
}
\examples{

library(RPANDA)

data(woodmouse)

alignment <- as.character(woodmouse) # nucleotidic alignment 

theta_estimator(alignment)

}

