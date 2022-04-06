\name{add.gts}
\alias{add.gts}

\title{
Geological time scale
}
\description{
Adds geological time scale (GTS) to plots.
}
\usage{
add.gts(thickness, quaternary = T, is.phylo = F,
        xpd.x = T, time.interval = 1, names = NULL, fill = T,
        cex = 1, padj = -0.5, direction = "rightwards")
}

\arguments{
  \item{thickness}{  
  numeric < 0. Define the thickness of the scale.
  }
  \item{quaternary}{
  bolean. Whether to merge Pleistocene and Holocene into Quaternary. Default is TRUE.
  }
  \item{is.phylo}{
  bolean. Whether the plot is a phylogeny or not. Default is FALSE.
  }
  \item{time.interval}{
  numeric. Define the minimum time interval (in million years) for the geological time scale. Default is 1 and displays ticks every million year but with numbers at every five million years.
  }
  \item{xpd.x}{
  bolean. Whether to expand the last period of the geological time scale before root age (mainly for tree). Default is TRUE.
  }
  \item{names}{
  a character vector with the names of geological periods (stages). Can be used to write abbreviations. Default is NULL and display full names (except for Quaternary and Pliocene). 
  }
  \item{fill}{
  bolean. If TRUE (default), backbground is alternatively filled with grey and white bands to distinguish geological periods. If FALSE, dashed lines are drawn to limit geological periods.
  }
  \item{cex}{
    numeric. Size of the names of geological periods.
  }
  \item{padj}{
    padj argument for the axis (see par() for more details).
  }
  \item{direction}{
  character. Direct the geological time scale. Can be either "rightwards" (default) of "leftwards" (NOT IMPLEMENTED YET). 
  }
}

\details{
  This function plots a geological times scale (GTS). It has been designed for adding GTS to plot of phylogeny, diversification rates and paleodiversity dynamics through time but can be used with any R plot. Time should be negative for other plots than phylogenies.
}

\value{
  Draws geological time scale on x axis.
}
\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep)
}

\author{Nathan Mazet}

\examples{
\dontrun{
# with a phylogeny

data("Cetacea")
# first plot to get the dimensions of the gts
plot(Cetacea, cex = 0.5, label.offset = 0.2, tip.color = "white")
add.gts(-11, quaternary = T, is.phylo = T, xpd.x = F,
        names = c("Q.", "Pli.", "Miocene", "Oligocene", "Eoc."))
# second plot to display the tree on the gts
par(new = T)
plot(Cetacea, cex = 0.5, label.offset = 0.2)
mtext("Time (Myrs)", side = 1, line = 3, at = 18)

# see Appendix S3 from Mazet et al. (in prep) for more examples.

}
}