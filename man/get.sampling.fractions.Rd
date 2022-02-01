\name{get.sampling.fractions}
\alias{get.sampling.fractions}

\title{
  Sampling fractions of subclades
}
\description{
  Provides the sampling fractions of a phylogenetic tree from a complete database.
}
\usage{
  get.sampling.fractions(phylo, data, clade.size = 5, plot = F,
                         lad = T, text.cex = 1, pch.cex = 0.8, ...)
}

\arguments{
  \item{phylo}{
    an object of type 'phylo' (see ape documentation)
  }
  \item{data}{
    an object of class "data.frame" containing a database of monophyletic groups for which potential shift can be tested. This database should be based on taxonomy, ecology or traits and must contain a column named "Species" with species names as in phylo. 
  }
  \item{clade.size}{
    numeric. Define the minimum number of species in a subgroup. Default is 5.
  }
  \item{plot}{
    bolean. If TRUE, the tree is plotted and testable nodes are highlighted with red dots. Default is FALSE.
  }
  \item{lad}{
    bolean. Define which way the tree should be represented if plot = T. If TRUE, the smallest clade is at the bottom plot. If FALSE, it is at the the top of the plot. Default is TRUE.
  }
    
  \item{text.cex}{
    numeric. Defines the size of the text in legend.
  }
  \item{pch.cex}{
   numeric. Defines the size of the red points at the crown of subclades.
  }
  \item{...}{
   
  }
}


\value{
  a dataframe with as many rows as nodes in the phylogeny with the following informations in columns:
  \item{nodes}{the node IDs}
  \item{data}{the name of the subclade from data}
  \item{f}{the sampling fraction for this subclade}
  \item{sp_in}{the number of species included in the tree}
  \item{sp_tt}{the number of species described in the data}
  \item{to_test}{the node IDs for nodes that are testable according to clade.size}
}
\details{
  All described species should be included to properly calculate sampling fractions. The example of Cetacea uses a taxonomic database but groups can be defined on geography or traits as soon as they are monophyletic. If the taxonomy of the studied group is difficult to establish (e.i. taxonomic uncertainty, etc.), a "fake" taxonomic database can be created with random species names (Gen1_sp1, Gen1_sp2, Gen2_sp1, etc.) to circumvent taxonomic difficulties. 

}
\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
}
\author{Nathan Mazet}
\seealso{
\code{\link{get.comb.shift}}, \code{\link{shift.estimates}}}

\examples{
# loading data

data("Cetacea")
data("taxo_cetacea")

# calculating sampling fractions
taxo_cetacea_no_genus <- taxo[names(taxo) != "Genus"]
f_cetacea <- get.sampling.fractions(phylo = Cetacea, lad = F,
                                    data = taxo_cetacea_no_genus,
                                    plot = T, cex = 0.3)

}