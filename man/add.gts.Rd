\name{add.gts}
\alias{add.gts}

\title{
Geological time scale
}
\description{
Adds geological time scale to plots.
}
\usage{
add.gts(Y1, quaternary = T, plot_dim = par("usr"), is.phylo = F,
              root.age = NULL, present = NULL, xpd.x = T,
              time.interval = 1, names = NULL, fill = T,
              cex = 1, padj = -0.5, direction = "rightwards")
}

\arguments{
  \item{thickness}{  
  numeric. Define the thickness of the scale.
  }
  \item{quaternary}{
  bolean. Whether to merge Pleistocene and Holocene into Quaternary. Default is TRUE.
  }
  \item{is.phylo}{
  bolean. Whether the plot is a plot of a phylogeny or not. Default is FALSE. If TRUE, present and root.age do not need to be mentionned.
  }
  \item{root.age}{
  numeric. Corresponds to the x coordinate of the more ancient time in the plot.
  }
  \item{present}{
  numeric. Corresponds to the x coordinate at present in the plot.
  }
  \item{xpd.x}{
  bolean. Whether to expand the last period of the geological time scale before root.age. Default is TRUE.
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
  This function plots a geological times scale (GTS). It has been designed for adding GTS to plot of phylogeny, diversification rates and paleodiversity dynamics through time but can be used with any R plot if present and root.age is correctly defined.
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
More to come soon...

# with a phylogeny

data("Cetacea")
plot(Cetacea, cex = 0.5, label.offset = 0.2, tip.color = "white")
add.gts(-11, quaternary = T, is.phylo = T, xpd.x = F,
        names = c("Q.", "Pli.", "Miocene", "Oligocene", "Eoc."))
par(new = T)
plot(Cetacea, cex = 0.5, label.offset = 0.2)
mtext("Time (Myrs)", side = 1, line = 3, at = 18)


# with rates 


# with diversity

}
}