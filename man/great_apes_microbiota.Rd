\name{great_apes_microbiota}
\alias{great_apes_microbiota}

\title{
Example of great apes microbiota
}
\description{
Phylogenetic tree of the great apes and example of three bacterial nucleotidic alignments (16S rRNA gene).
}
\usage{
data(great_apes_microbiota)
}

\author{Benoit Perez-Lamarque}
\references{
Ochman H, Worobey M, Kuo CH, Ndjango JBN, Peeters M, Hahn BH, Hugenholtz P (2010). Evolutionary relationships of wild hominids recapitulated by gut microbial communities. PLoS Biology, 8(11), 3-10. https://doi.org/10.1371/journal.pbio.1000546.

Sanders JG, Powell S, Kronaue DJC, Vasconcelos HL, Fredrickson ME, Pierce NE (2014). Data from: Stability and phylogenetic correlation in gut microbiota: lessons from ants and apes. Dryad Digital Repository, https://doi.org/10.5061/dryad.023s6.

Perez-Lamarque B, Morlon H (2019). Characterizing symbiont inheritance during host-microbiota evolution: Application to the great apes gut microbiota. Molecular Ecology Resources 19:1659-1671.
	}

\seealso{\code{\link{example_great_apes_microbiota}},
\code{\link{HOME_model}}
}


\examples{
\dontshow{.old_wd <- setwd(tempdir())}
data(great_apes_microbiota, envir = environment())
\dontshow{setwd(.old_wd)}
}
