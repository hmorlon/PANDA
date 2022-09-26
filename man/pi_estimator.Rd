\name{pi_estimator}
\alias{pi_estimator}

\title{
Compute nucleotidic diversity (Pi estimator)
}

\description{
This function computes the Pi estimator of genetic diversity (Nei and Li, 1979) while controlling for the presence of gaps in the alignment (Ferretti et al, 2012), frequent in barcoding datasets. 
}

\usage{
pi_estimator(sequences)
}

\arguments{

\item{sequences}{
a matrix representing the nucleotidic alignment of all the sequences present in the phylogenetic tree. 
}

}


\value{
An estimate of genetic diversity 
}

\references{
Nei M & Li WH, Mathematical model for studying genetic variation in terms of restriction endonucleases, 1979, Proc. Natl. Acad. Sci. USA.

Ferretti L, Raineri E, Ramos-Onsins S. 2012. Neutrality tests for sequences with missing data. Genetics 191: 1397–1401.

Perez-Lamarque B, Öpik M, Maliet O, Silva A, Selosse M-A, Martos F, and Morlon H. 2022. Analysing diversification dynamics using barcoding data: The case of an obligate mycorrhizal symbiont, Molecular Ecology, 31:3496–512.
}
\author{
Ana C. Afonso Silva & Benoît Perez-Lamarque}
\seealso{
\code{\link{theta_estimator}}
\code{\link{delineate_phylotypes}}
}
\examples{

library(RPANDA)

data(woodmouse)

alignment <- as.character(woodmouse) # nucleotidic alignment 

pi_estimator(alignment)

}

