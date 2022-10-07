\name{plot.phylo.comb}
\alias{plot.phylo.comb}

\title{
  Plot shifts of diversifcation on a phylogeny
}
\description{
  Plots the phylogeny with colored branches according to shifts of diversification.
}

\usage{
plot.phylo.comb(phylo, data, sampling.fractions, shift.res = NULL,
                combi, backbone.option = "crown.shift",
                main = NULL, col.sub = NULL, col.bck = "black",
                lty.bck = 1, tested_nodes = F, lad = T,
                leg = T, text.cex = 1, pch.cex = 1, \dots)
}
\arguments{
  \item{phylo}{
    an object of type 'phylo' (see ape documentation)
  }
  \item{data}{
    a data.frame containing a database of monophyletic groups for which potential shifts can be tested. This database should be based on taxonomy, ecology or traits and must contain a column named "Species" with species names as in phylo. 
  }
  \item{sampling.fractions}{
    the output resulting from get.sampling.fractions.
  }
  \item{shift.res}{
    the output resulting from shift.estimates or NULL (default). This latter case allows to represent combinations only from the output of \code{\link{get.comb.shift}} by specifying the combination (see argument combi).
  }
  \item{combi}{
    character or numeric. If shift.res is provided, this argument is a numeric and corresponds to the rank of the combination in the global comparison (shift.res$total). If shift.res is NULL, this argument should be a character giving a combination of node IDs as in get.comb.shift output. This latter manner to specify combination allows to visualize a combination of shifts before having results.
  }
  \item{backbone.option}{
    type of the backbone analysis:
    \itemize{ 
      \item"stem.shift": the stems of subclades are included in subclade analyses;
      \item"crown.shift": the stems of subclades are included in the backbone analysis (Default).
    }
  }
  \item{main}{
    Character. The name of the plot. Default is NULL and the combination rank with AICc will be printed if shift.res is not NULL.
  }
  \item{col.sub}{
    character. A vector to specify colors of subclade(s). Can be let NULL (see details).
  }
  \item{col.bck}{
    character. A vector to specify colors of backbone(s). Default is "black" for simple backbone (see details).
  }
  \item{lad}{
    bolean. Allows to ladderize the tree.
  }
  \item{leg}{
    bolean. If TRUE, legend of the selected combination is added to the plot with names from data and best model names. Default is TRUE. The position is automatically adjusted in function of lad argument.
  }
  \item{lty.bck}{
    numeric. Define lty for the backbone.
  }

  \item{tested_nodes}{
    bolean. If TRUE, all the tested nodes are highlighted by a red point.
  }
  \item{text.cex}{
    numeric. Define the size of legend text.
  }
  \item{pch.cex}{
    numeric. Define the size of points if tested_nodes = TRUE
  }
  \item{...}{
    further arguments to be passed to plot or to plot.phylo.
  }
}

\value{
plot the phylogeny and returns the same invisible object as plot.phylo.
}
\details{
  If col.sub is not specified, color vector for subclades is c(c(brewer.pal(8, "Dark2"),brewer.pal(8, "Set1"),"darkmagenta","dodgerblue2", "orange", "forestgreen")). For multiple backbone, default vector is c("blue4", "orange4", "red4", "grey40", "coral4", "deeppink4", "khaki4", "darkolivegreen", "darkslategray",”black”).
  ... allows to set different graphical parameters from plot.phylo such as cex for size of tip labels or edge.width for the thickness of the phylogeny edges.

}
\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
}
\author{Nathan Mazet}
\seealso{
\code{\link{shift.estimates}}}

\examples{
\donttest{
# loading data
data("Cetacea")
data("taxo_cetacea")

taxo_cetacea_no_genus <- taxo_cetacea[names(taxo_cetacea) != "Genus"]

# main procedure
f_cetacea <- get.sampling.fractions(phylo = Cetacea, lad = FALSE,
                                    data = taxo_cetacea_no_genus,
                                    plot = TRUE, cex = 0.3)

comb.shift_cetacea <- get.comb.shift(phylo = Cetacea,
                                     data = taxo_cetacea_no_genus,
                                     sampling.fractions = f_cetacea,
                                     Ncores = 4)

# use of plot.phylo.comb
# without shift.estimates results but with comb.shift_cetacea

plot.phylo.comb(phylo = Cetacea,
                data = taxo_cetacea,
                sampling.fractions = f_cetacea,
                combi = comb.shift_cetacea[15],
                label.offset = 0.3,
                main = "", lad = FALSE ,cex = 0.4)
}
}