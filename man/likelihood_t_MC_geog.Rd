\name{likelihood_t_MC_geog}
\alias{likelihood_t_MC_geog}

\title{
Likelihood of a dataset under the matching competition model with biogeography.
}
\description{
Computes the likelihood of a dataset under the matching competition model with specified \code{sigma2} and \code{S} values and with a \code{geography.object} formed using \code{\link{CreateGeoObject}}. 
}
\usage{
likelihood_t_MC_geog(phylo, data, par,geography.object)
}

\arguments{
  \item{phylo}{  
an object of type 'phylo' (see ape documentation)
}
  
  \item{data}{
a named vector of continuous data with names corresponding to \code{phylo$tip.label}
}
 \item{par}{
a vector listing a value for \code{sig2} and \code{S} (parameters of the matching competition model), in that order
}
  \item{geography.object}{
a list specifying sympatry (= 1) and allopatry (= 0) for each extant lineages at each internode interval (see details)}

}

\details{
	When specifying \code{par}, \code{sig2} must be listed before \code{S}.
	
	The \code{geography.object} should be a list of matrices with row and column
names that match the \code{$geography.object} object created with the \code{\link{CreateGeoObject}} function.

}

\note{
	S must be negative (if it is positive, the likelihood function will multiply input by -1).
}

\value{
  the negative log-likelihood value of the dataset (accordingly, the negative of the output should be recorded as the likelihood), given the phylogeny, \code{sig2} and \code{S} values, and \code{geography.object}.
}
\references{
	Drury, J., Clavel, J., and Morlon, H. (submitted) Estimating the effect of competition on trait evolution using maximum likelihood inference.

	Nuismer, S. & Harmon, L. 2015. Predicting rates of interspecific interaction from phylogenetic trees. \emph{Ecology Letters} 18:17-27.
}

\author{Jonathan Drury \href{mailto:drury@biologie.ens.fr}{drury@biologie.ens.fr}

		Julien Clavel
}

\seealso{
\code{\link{fit_t_comp}}
\code{\link{CreateGeoObject}}
\code{\link{likelihood_t_MC}}

}


\examples{
data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data
geography.object <-  Anolis.data$geography.object

# Compute the likelihood with geography using ML parameters for fit without geography
par <- c(0.0003139751, -0.06387258)
lh <- -likelihood_t_MC_geog(phylo,pPC1,par,geography.object)
}
